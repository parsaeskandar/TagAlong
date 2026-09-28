// Coordinate translation using ONLY a GBZ and a gbwt::SequenceLocate (.sri).
//
// The existing implementation (coordinate_translation.cpp + pangenome_server.cpp)
// needs five structures: an 81 GB RLBWT r-index, a 20 GB sampled tag array, an
// 11 GB GBWT FastLocate, Table 1 and Table 2 (~120 GB, ~11 min to load). All of
// that scaffolding exists for ONE reason: the RLBWT was the only way to turn a
// graph node into a base offset on a haplotype. That forced the anchor search
// (first/last common node, a uniqueness pass, a colinearity gate), the backward
// extended search, and a second pass against the reverse sequence to catch
// inversions.
//
// SequenceLocate answers that question directly, in both directions:
//   locateForward(state)      node  -> (sequence, base offset) for every visit
//   locateSequence(seq, bp)   base  -> (node, offset) + the GBWT position
//
// So this file needs no anchors, no tag array, no tables:
//
//   1. metadata.findFragment()     resolve the source contig interval (was Table 1)
//   2. locateSequence() ONCE       enter the graph at the interval start
//   3. LF-walk the source path     one node at a time, accumulating base offsets
//   4. locateForward() per node    every haplotype at that node, with its base offset
//   5. emit one block per node     strand from comparing Path::is_reverse
//
// Inversions need no special case: querying the FORWARD node returns each visit
// exactly once, a forward traversal as the forward sequence and a backward one as
// the reverse sequence, so relative strand is a comparison of Path::is_reverse.
//
// Built separately from liftover_ext because it links a DIFFERENT gbwt (the one
// with SequenceLocate) than the rest of the project; mixing them in one binary
// would be an ODR violation.

#include <gbwt/sequence_locate.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace seqtrans {

/// One node's worth of correspondence: a source range and the target range it
/// maps to. Blocks are the natural output here -- every node gives a whole
/// range at once, where the old code produced one pair per base and folded
/// afterwards.
struct Block {
    std::string target_contig;
    size_t src_start = 0;      ///< contig coordinates, half-open
    size_t src_end   = 0;
    size_t tgt_start = 0;
    size_t tgt_end   = 0;
    char   strand    = '+';
    size_t tgt_path_id = 0;
};

struct Stats {
    size_t fragments        = 0;
    size_t nodes_walked     = 0;
    size_t locate_sequence  = 0;   ///< should be ~1 per fragment
    size_t locate_forward   = 0;   ///< one per node
    size_t visits_enumerated= 0;   ///< total locateNext cost
    size_t blocks           = 0;
    size_t self_check_fail  = 0;   ///< walked offset != index's own answer
    size_t target_walks     = 0;   ///< WALK mode: target sequences walked
    size_t target_walk_steps= 0;   ///< WALK mode: LF steps along the targets
    size_t verify_checked   = 0;   ///< nodes re-checked against locateForward
    size_t verify_mismatch  = 0;   ///< ... that disagreed
    double ms               = 0.0;
};

/// SCAN calls locateForward on every source node: simple, and every node's
/// answer comes straight from the index. WALK calls it only until it finds the
/// target, then follows that target with LF and joins on node identity -- the
/// same trick that took surjection from 310 s to 0.37 s. VERIFY runs WALK and
/// re-checks a sample of nodes against SCAN.
enum class Mode { SCAN, WALK, VERIFY };

class Translator {
public:
    Translator(const gbwtgraph::GBZ& gbz, const gbwt::SequenceLocate& sri)
        : gbz_(gbz), sri_(sri) {}

    /// Path ids whose (sample, haplotype) match `name`, which may be
    /// "HG02257#1" or a full contig "GRCh38#0#chr8".
    std::vector<gbwt::size_type> resolve(const std::string& name) const {
        std::vector<gbwt::size_type> out;
        const gbwt::Metadata& md = gbz_.index.metadata;
        for (gbwt::size_type pid = 0; pid < md.paths(); pid++) {
            if (path_matches(pid, name)) out.push_back(pid);
        }
        return out;
    }

    /// Contig name as the rest of the stack spells it: sample#hap#contig.
    std::string contig_name(gbwt::size_type path_id) const {
        gbwt::FullPathName fp = gbz_.index.metadata.fullPath(path_id);
        return fp.sample_name + "#" + std::to_string(fp.haplotype) + "#" + fp.contig_name;
    }

    /// Translate [start, end) on `src_name` to every path in `tgt_ids`.
    std::vector<Block> translate(const std::string& src_name,
                                 size_t start, size_t end,
                                 const std::vector<gbwt::size_type>& tgt_ids,
                                 Stats& st, Mode mode = Mode::WALK) const {
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<Block> out;
        if (end <= start) return out;

        std::unordered_set<gbwt::size_type> tgt_set(tgt_ids.begin(), tgt_ids.end());

        // ── 1. source fragments (replaces Table 1) ─────────────────────────
        // A contig is split into several GBWT paths; walk the ones the interval
        // touches, in contig order.
        for (gbwt::size_type src_pid : fragments_for(src_name, start, end)) {
            gbwt::FullPathName fp = gbz_.index.metadata.fullPath(src_pid);
            const gbwt::size_type src_seq = gbwt::Path::encode(src_pid, false);
            const size_t frag_len = sri_.sequenceLength(src_seq);
            const size_t frag_off = fp.offset;

            // Clip the request to this fragment, in fragment-local coordinates.
            const size_t lo = (start > frag_off) ? (start - frag_off) : 0;
            const size_t hi = std::min(end - frag_off, frag_len);
            if (lo >= hi) continue;
            st.fragments++;

            if (mode != Mode::SCAN) {
                walk_fragment(out, st, src_pid, src_seq, frag_off, lo, hi,
                              tgt_set, mode == Mode::VERIFY);
                continue;
            }

            // ── 2. enter the graph ONCE ────────────────────────────────────
            // locateSequence is the expensive query (v2.0: mean 673us, worst
            // 161ms at the far end of an 18 Mb gap), so it runs once per
            // fragment and the rest of the interval is reached by LF.
            gbwt::edge_type gpos;
            gbwt::edge_type base = sri_.locateSequence(src_seq, lo, gpos);
            st.locate_sequence++;
            if (base == gbwt::invalid_edge() || gpos == gbwt::invalid_edge()) continue;

            // Start of the node containing `lo`, in fragment coordinates.
            size_t node_start = lo - base.second;

            // ── 3. walk the source path one node at a time ─────────────────
            while (gpos.first != gbwt::ENDMARKER && node_start < hi) {
                const gbwt::node_type node = gpos.first;
                const size_t L = node_len(node);
                st.nodes_walked++;

                // ── 4. every haplotype at this node, WITH base offsets ─────
                // Query the node in the SOURCE's traversal orientation, not the
                // graph-forward one. Every hit then traverses this node the same
                // oriented way the source does, so:
                //   - the source's own entry is its FORWARD sequence, and its
                //     offset is already a path offset (no frame conversion);
                //   - a target appearing as its REVERSE sequence is one that
                //     runs through here backwards, i.e. an inversion;
                //   - LF on any hit advances in the same direction as the
                //     source, which is what keeps the WALK below parallel.
                // Each visit still appears exactly once, so nothing is double
                // counted. (The note's "query the forward node" recipe is for
                // starting from a graph position, not from a source traversal.)
                gbwt::size_type first = gbwt::SequenceLocate::NO_POSITION;
                gbwt::SearchState state = sri_.find(node, first);
                st.locate_forward++;
                if (state.empty()) { advance(gpos, node_start, L); continue; }

                std::vector<std::pair<gbwt::size_type, gbwt::size_type>> hits =
                    sri_.locateForward(state, first);
                st.visits_enumerated += hits.size();

                // The source is in its own result; use it for the strand
                // reference and as a free check that the walk and the index
                // agree on where we are.
                bool have_src = false;
                for (const auto& h : hits) {
                    if (h.first != src_seq) continue;        // our forward sequence
                    if (h.second != node_start) continue;    // a different visit
                    have_src = true;
                    break;
                }
                if (!have_src) st.self_check_fail++;

                // ── 5. one block per (node, target visit) ──────────────────
                for (const auto& h : hits) {
                    const gbwt::size_type seq = h.first;
                    const gbwt::size_type pid = gbwt::Path::id(seq);
                    if (!tgt_set.count(pid)) continue;
                    if (pid == src_pid) continue;            // identity
                    // Queried in the source's orientation, so a reverse
                    // sequence here IS an inverted correspondence.
                    emit(out, st, pid, seq, h.second, L, node_start,
                         frag_off, lo, hi, gbwt::Path::is_reverse(seq));
                }
                advance(gpos, node_start, L);
            }
        }
        st.ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        return out;
    }

    /// Merge adjacent blocks that keep a constant source->target offset, so the
    /// output is a handful of alignment blocks rather than one entry per node.
    static std::vector<Block> fold(std::vector<Block> in) {
        if (in.empty()) return in;
        std::sort(in.begin(), in.end(), [](const Block& a, const Block& b) {
            if (a.tgt_path_id != b.tgt_path_id) return a.tgt_path_id < b.tgt_path_id;
            if (a.strand != b.strand) return a.strand < b.strand;
            return a.src_start < b.src_start;
        });
        std::vector<Block> out;
        for (const Block& b : in) {
            if (!out.empty()) {
                Block& p = out.back();
                const bool same = p.tgt_path_id == b.tgt_path_id && p.strand == b.strand;
                if (same && p.src_end == b.src_start &&
                    ((b.strand == '+' && p.tgt_end == b.tgt_start) ||
                     (b.strand == '-' && b.tgt_end == p.tgt_start))) {
                    p.src_end = b.src_end;
                    if (b.strand == '+') p.tgt_end = b.tgt_end;
                    else                 p.tgt_start = b.tgt_start;
                    continue;
                }
            }
            out.push_back(b);
        }
        return out;
    }

private:
    const gbwtgraph::GBZ& gbz_;
    const gbwt::SequenceLocate& sri_;

    /// Emit one node's correspondence, CLIPPED to the requested interval.
    /// A node is a whole unit in the index but the request is arbitrary, so the
    /// first and last node of an interval usually stick out at the ends. The
    /// clip has to move the target edge by the same amount -- and from the
    /// OPPOSITE end when the block is inverted, since source and target run in
    /// opposite directions there.
    void emit(std::vector<Block>& out, Stats& st,
              gbwt::size_type pid, gbwt::size_type seq, size_t fwd_start,
              size_t L, size_t node_start, size_t frag_off,
              size_t lo, size_t hi, bool inverted) const {
        size_t s0 = node_start, s1 = node_start + L;          // fragment coords
        const size_t cut_front = (lo > s0) ? (lo - s0) : 0;
        const size_t cut_back  = (s1 > hi) ? (s1 - hi) : 0;
        if (cut_front + cut_back >= L) return;                // nothing left
        s0 += cut_front; s1 -= cut_back;

        const bool rev  = gbwt::Path::is_reverse(seq);
        const size_t tl = sri_.sequenceLength(seq);
        // locateForward reports the node's start in the sequence's own
        // orientation; put it back on the path.
        size_t t0 = rev ? (tl - fwd_start - L) : fwd_start;
        size_t t1 = t0 + L;
        if (inverted) { t0 += cut_back;  t1 -= cut_front; }   // mirrored
        else          { t0 += cut_front; t1 -= cut_back;  }

        gbwt::FullPathName tfp = gbz_.index.metadata.fullPath(pid);
        Block b;
        b.target_contig = contig_name(pid);
        b.tgt_path_id   = pid;
        b.src_start = frag_off + s0;  b.src_end = frag_off + s1;
        b.tgt_start = tfp.offset + t0; b.tgt_end = tfp.offset + t1;
        b.strand    = inverted ? '-' : '+';
        out.push_back(b);
        st.blocks++;
    }

    size_t node_len(gbwt::node_type node) const {
        return gbz_.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(node));
    }

    void advance(gbwt::edge_type& gpos, size_t& node_start, size_t L) const {
        node_start += L;
        gpos = gbz_.index.LF(gpos);     // navigate with the GBWT position, never
    }                                    // with the (node, offset_in_node) value

    /// WALK mode for one source fragment.
    ///
    /// SCAN asks the index about every source node, which costs one locateNext
    /// per VISIT -- ~146 visits per node on HPRC v2.0, and we keep the one or
    /// two that belong to the target. Instead: walk the source once, ask
    /// locateForward only until a target turns up, then follow that target with
    /// LF and join the two walks on node identity. Both walks accumulate base
    /// offsets by summing node lengths, so no further index queries are needed.
    ///
    /// Divergence needs no special handling: a node only the source has simply
    /// finds no partner, and an insertion on the target is skipped by the join.
    void walk_fragment(std::vector<Block>& out, Stats& st,
                       gbwt::size_type src_pid, gbwt::size_type src_seq,
                       size_t frag_off, size_t lo, size_t hi,
                       const std::unordered_set<gbwt::size_type>& tgt_set,
                       bool verify) const {
        gbwt::edge_type gpos;
        gbwt::edge_type base = sri_.locateSequence(src_seq, lo, gpos);
        st.locate_sequence++;
        if (base == gbwt::invalid_edge() || gpos == gbwt::invalid_edge()) return;

        // ── source walk: (node, fragment-local start), in path order ───────
        std::vector<std::pair<gbwt::node_type, size_t>> src_nodes;
        size_t node_start = lo - base.second;
        while (gpos.first != gbwt::ENDMARKER && node_start < hi) {
            const size_t L = node_len(gpos.first);
            src_nodes.emplace_back(gpos.first, node_start);
            st.nodes_walked++;
            node_start += L;
            gpos = gbz_.index.LF(gpos);
        }
        if (src_nodes.empty()) return;

        // ── find the target(s): locateForward until one shows up ───────────
        // Every target visit at that node starts its own walk, so a region
        // present in several copies is followed in all of them.
        struct TWalk { gbwt::size_type seq; gbwt::edge_type pos; size_t bp; };
        std::vector<TWalk> walks;
        size_t probe = 0;
        for (; probe < src_nodes.size() && walks.empty(); probe++) {
            // Probe in the SOURCE's orientation. Every hit then traverses
            // this node the way the source does, so LF advances all of them in
            // the same direction -- including a target running through the
            // region backwards, which arrives as its reverse sequence.
            // Anchoring on the graph-forward node sent anti-parallel targets
            // off the wrong way: 21,393 LF steps, 1 node matched.
            const gbwt::node_type snode = src_nodes[probe].first;
            gbwt::size_type first = gbwt::SequenceLocate::NO_POSITION;
            gbwt::SearchState state = sri_.find(snode, first);
            st.locate_forward++;
            if (state.empty()) continue;
            auto hits = sri_.locateForward(state, first);
            st.visits_enumerated += hits.size();
            for (size_t i = 0; i < hits.size(); i++) {
                const gbwt::size_type seq = hits[i].first;
                if (!tgt_set.count(gbwt::Path::id(seq))) continue;
                if (gbwt::Path::id(seq) == src_pid) continue;
                // entry i of the state's range IS the GBWT record offset
                walks.push_back({seq,
                                 gbwt::edge_type(snode, state.range.first + i),
                                 hits[i].second});
                st.target_walks++;
            }
        }
        if (walks.empty()) return;

        // ── walk each target, recording where it meets each node ───────────
        // node -> (forward start in that sequence). Budget is a safety valve;
        // the walk normally ends when it runs past the source's span.
        const size_t budget = std::max<size_t>(1u << 16, src_nodes.size() * 8);
        std::unordered_map<gbwt::node_type, std::vector<std::pair<size_t,size_t>>> tmap;
        for (size_t w = 0; w < walks.size(); w++) {
            gbwt::edge_type tp = walks[w].pos;
            size_t tbp = walks[w].bp;
            size_t steps = 0;
            while (tp.first != gbwt::ENDMARKER && steps < budget) {
                tmap[tp.first].emplace_back(w, tbp);
                tbp += node_len(tp.first);
                tp = gbz_.index.LF(tp);
                steps++; st.target_walk_steps++;
                if (steps > src_nodes.size() * 4) break;   // well past the source
            }
        }

        // ── join on node identity ──────────────────────────────────────────
        for (size_t k = 0; k < src_nodes.size(); k++) {
            const gbwt::node_type snode = src_nodes[k].first;
            const size_t L = node_len(snode);
            // Two independent things decide the strand, and BOTH matter:
            //
            //   flip  -- the target walk met this node in the opposite
            //            orientation, i.e. a locally inverted segment inside an
            //            otherwise parallel alignment;
            //   is_reverse(seq) -- the walk is following the target's REVERSE
            //            sequence, so its path coordinate runs backwards
            //            relative to the walk (a target crossing the whole
            //            region the other way round).
            //
            // Either one alone inverts the correspondence; both together cancel
            // out. Hence XOR. Dropping the flip loses the inverted core of an
            // otherwise collinear alignment (the boundary fixture); ignoring
            // is_reverse mislabels an anti-parallel target as '+'.
            for (int flip = 0; flip < 2; flip++) {
                const gbwt::node_type key = flip
                    ? gbwt::Node::encode(gbwt::Node::id(snode),
                                         !gbwt::Node::is_reverse(snode))
                    : snode;
                auto it = tmap.find(key);
                if (it == tmap.end()) continue;
                for (const auto& hit : it->second) {
                    const TWalk& tw = walks[hit.first];
                    const bool inverted =
                        ((flip != 0) != gbwt::Path::is_reverse(tw.seq));
                    // hit.second is the node's start in that sequence's own
                    // orientation, which is what emit() expects.
                    emit(out, st, gbwt::Path::id(tw.seq), tw.seq, hit.second,
                         L, src_nodes[k].second, frag_off, lo, hi, inverted);
                }
            }
        }

        // ── optional spot-check against SCAN ───────────────────────────────
        if (verify) {
            for (size_t k = 0; k < src_nodes.size(); k += 512) {
                gbwt::size_type first = gbwt::SequenceLocate::NO_POSITION;
                gbwt::SearchState state = sri_.find(src_nodes[k].first, first);
                if (state.empty()) continue;
                auto hits = sri_.locateForward(state, first);
                size_t expect = 0;
                for (const auto& h : hits)
                    if (tgt_set.count(gbwt::Path::id(h.first)) &&
                        gbwt::Path::id(h.first) != src_pid) expect++;
                size_t got = 0;
                for (int flip = 0; flip < 2; flip++) {
                    const gbwt::node_type key = flip
                        ? gbwt::Node::encode(gbwt::Node::id(src_nodes[k].first),
                                             !gbwt::Node::is_reverse(src_nodes[k].first))
                        : src_nodes[k].first;
                    auto it2 = tmap.find(key);
                    if (it2 != tmap.end()) got += it2->second.size();
                }
                st.verify_checked++;
                if (got != expect) st.verify_mismatch++;
            }
        }
    }

    bool path_matches(gbwt::size_type pid, const std::string& name) const {
        gbwt::FullPathName fp = gbz_.index.metadata.fullPath(pid);
        const std::string full = fp.sample_name + "#" + std::to_string(fp.haplotype)
                               + "#" + fp.contig_name;
        if (full == name) return true;
        const std::string hap = fp.sample_name + "#" + std::to_string(fp.haplotype);
        return hap == name;
    }

    /// Fragments of `src_name` overlapping [start, end), in contig order.
    std::vector<gbwt::size_type> fragments_for(const std::string& src_name,
                                               size_t start, size_t end) const {
        std::vector<std::pair<size_t, gbwt::size_type>> hits;
        const gbwt::Metadata& md = gbz_.index.metadata;
        for (gbwt::size_type pid = 0; pid < md.paths(); pid++) {
            if (!path_matches(pid, src_name)) continue;
            gbwt::FullPathName fp = md.fullPath(pid);
            const size_t flen = sri_.sequenceLength(gbwt::Path::encode(pid, false));
            if (fp.offset < end && start < fp.offset + flen) hits.emplace_back(fp.offset, pid);
        }
        std::sort(hits.begin(), hits.end());
        std::vector<gbwt::size_type> out;
        out.reserve(hits.size());
        for (const auto& h : hits) out.push_back(h.second);
        return out;
    }
};

} // namespace seqtrans
