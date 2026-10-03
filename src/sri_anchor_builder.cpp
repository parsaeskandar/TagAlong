// Surjection anchor construction from a SequenceLocate (.sri) index alone.
//
// The existing builder needs FastLocate (which visits exist, and their record
// offsets) plus the RLBWT tag array (base offsets), then zips the two sorted
// lists positionally -- a pairing that silently misaligns when the indexes
// disagree, which on HPRC v2.0 they do. SequenceLocate stores offsets in base
// pairs, so one decompressSA() answers both halves at once: sa[k] packs
// (sequence id, bp offset), and since gbwt::find returns an inclusive range
// based at 0, the index k IS the GBWT record offset. No zip, no LF walk.

#include <gbwt/gbwt.h>
#include <gbwt/sequence_locate.h>
#include <gbwtgraph/gbz.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace srianchor {

/// The subset of panindexer::SourceMapping the anchor search actually reads.
struct SourceMapping {
    int64_t node_id = 0;
    bool    is_reverse = false;
    size_t  read_begin_offset = 0;
    size_t  read_end_offset = 0;
};

/// One anchor, carrying the fields the A/B signature compares.
struct Anchor {
    size_t read_begin = 0;
    size_t read_end = 0;
    size_t path_offset = 0;          ///< forward base offset on the target path
    gbwt::node_type edge_node = 0;   ///< gbwt::edge_type.first
    size_t edge_offset = 0;          ///< gbwt::edge_type.second (offset in record)
    size_t source_mapping = 0;
};

struct Stats {
    size_t mappings = 0;
    size_t distinct_nodes = 0;
    size_t sa_calls = 0;
    size_t sa_entries = 0;
    size_t anchors = 0;
    size_t touched_paths = 0;
    bool   reverse_strand = false;
    double ms = 0.0;
};

/// One target's result from build_multi(): exactly what build() returns for
/// that target alone.
struct TargetResult {
    std::vector<Anchor> anchors;
    size_t best_pid = 0;
    bool reverse_strand = false;
    size_t touched_paths = 0;
};

class Builder {
public:
    Builder(const gbwtgraph::GBZ& gbz, const gbwt::SequenceLocate& sri)
        : gbz_(gbz), sri_(sri) {}

    /// Resolve a full path name or a haplotype prefix to GBWT path ids, using
    /// the same exact-then-prefix rule as the TranslationTable1 lookup.
    std::vector<size_t> resolve(const std::string& name) const {
        std::vector<size_t> exact, by_hap;
        const gbwt::Metadata& md = gbz_.index.metadata;
        for (size_t pid = 0; pid < md.paths(); pid++) {
            gbwt::FullPathName fp = md.fullPath(pid);
            const std::string hap = fp.sample_name + "#" + std::to_string(fp.haplotype);
            if (hap + "#" + fp.contig_name == name) exact.push_back(pid);
            else if (hap == name)                   by_hap.push_back(pid);
        }
        return exact.empty() ? by_hap : exact;
    }

    /// resolve() for many names in ONE pass over the metadata; resolving each
    /// separately walks all ~95k paths per name.
    std::vector<std::vector<size_t>> resolve_many(const std::vector<std::string>& names) const {
        std::unordered_map<std::string, std::vector<size_t>> where;
        for (size_t i = 0; i < names.size(); i++) where[names[i]].push_back(i);
        std::vector<std::vector<size_t>> exact(names.size()), by_hap(names.size());
        const gbwt::Metadata& md = gbz_.index.metadata;
        for (size_t pid = 0; pid < md.paths(); pid++) {
            gbwt::FullPathName fp = md.fullPath(pid);
            const std::string hap = fp.sample_name + "#" + std::to_string(fp.haplotype);
            auto e = where.find(hap + "#" + fp.contig_name);
            if (e != where.end()) for (size_t i : e->second) exact[i].push_back(pid);
            auto h = where.find(hap);
            if (h != where.end()) for (size_t i : h->second) by_hap[i].push_back(pid);
        }
        for (size_t i = 0; i < names.size(); i++) {
            if (exact[i].empty()) exact[i] = std::move(by_hap[i]);
        }
        return exact;
    }

    std::string path_name(size_t pid) const {
        gbwt::FullPathName fp = gbz_.index.metadata.fullPath(pid);
        return fp.sample_name + "#" + std::to_string(fp.haplotype) + "#" + fp.contig_name;
    }

    /// Build anchors for the best-matching target path among `target_pids`.
    ///
    /// Only sequence id 2*pid is matched. That is complete because a read on the
    /// target's reverse strand traverses every node flipped, so flipping the
    /// read's own mappings makes it match the target's forward sequence -- which
    /// is also why forwardStart() needs no mirroring here.
    std::vector<Anchor> build(const std::vector<SourceMapping>& mappings,
                              const std::vector<size_t>& target_pids,
                              Stats& st, size_t* best_pid_out = nullptr) const {
        const auto t0 = std::chrono::steady_clock::now();
        st = Stats{};
        st.mappings = mappings.size();
        if (mappings.empty() || target_pids.empty()) return {};

        std::unordered_set<size_t> want(target_pids.begin(), target_pids.end());

        // Probe the strand before scanning. Scanning forward first and only
        // then flipping costs a wasted full pass on every reverse-strand read;
        // the probe settles it in one or two decompressSA calls.
        const int strand = detect_strand(mappings, want, st);
        if (strand < 0) {
            st.ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
            return {};
        }
        st.reverse_strand = (strand == 1);

        std::unordered_map<size_t, std::vector<Anchor>> by_path;
        if (strand == 0) {
            scan(mappings, want, by_path, st);
        } else {
            scan(reverse_strand_mappings(mappings), want, by_path, st);
        }

        st.touched_paths = by_path.size();
        TargetResult r;
        finalize(by_path, r);
        if (best_pid_out) *best_pid_out = r.best_pid;
        st.anchors = r.anchors.size();
        st.ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        return std::move(r.anchors);
    }

    /// build() for many targets in one pass. `targets[t]` is target t's path-id
    /// set (as resolve() returns); result t is exactly build(mappings, targets[t]).
    ///
    /// decompressSA enumerates every haplotype at a node however many we keep,
    /// so one pass costs the same for 1 target or 50. The read-orientation pass
    /// serves strand detection and every forward-strand target; the flipped pass
    /// runs once, shared, only if some target is on the reverse strand.
    std::vector<TargetResult> build_multi(const std::vector<SourceMapping>& mappings,
                                          const std::vector<std::vector<size_t>>& targets,
                                          Stats& st) const {
        const auto t0 = std::chrono::steady_clock::now();
        st = Stats{};
        st.mappings = mappings.size();
        const size_t T = targets.size();
        std::vector<TargetResult> res(T);
        if (mappings.empty() || T == 0) return res;

        // pid -> the targets it belongs to (a target may be requested twice).
        std::unordered_map<size_t, std::vector<uint32_t>> owners;
        for (size_t t = 0; t < T; t++)
            for (size_t pid : targets[t]) owners[pid].push_back(static_cast<uint32_t>(t));
        if (owners.empty()) return res;

        // -1 undecided, 0 forward, 1 reverse -- per target, as detect_strand.
        std::vector<int> strand(T, -1);
        size_t undecided = 0;
        for (size_t t = 0; t < T; t++) if (!targets[t].empty()) undecided++;
        std::vector<std::unordered_map<size_t, std::vector<Anchor>>> by_path(T);

        // ── read-orientation pass ────────────────────────────────────────
        struct Entry { size_t pid; bool rev; size_t path_offset; size_t k; };
        std::unordered_map<gbwt::node_type, std::vector<Entry>> cache;
        std::unordered_map<uint32_t, bool> node_has_fwd;   // per node, per target
        for (size_t i = 0; i < mappings.size(); i++) {
            const gbwt::node_type node = gbwt::Node::encode(
                static_cast<gbwt::size_type>(mappings[i].node_id), mappings[i].is_reverse);
            auto it = cache.find(node);
            if (it == cache.end()) {
                std::vector<Entry> entries;
                std::vector<gbwt::size_type> sa = sri_.decompressSA(node);
                st.sa_calls++; st.sa_entries += sa.size(); st.distinct_nodes++;
                for (size_t k = 0; k < sa.size(); k++) {
                    const gbwt::size_type seq = sri_.seqId(sa[k]);
                    const size_t pid = gbwt::Path::id(seq);
                    if (!owners.count(pid)) continue;
                    const bool rev = gbwt::Path::is_reverse(seq);
                    entries.push_back({pid, rev,
                        rev ? 0 : sri_.forwardStart(seq, sri_.seqOffset(sa[k]), node), k});
                }
                it = cache.emplace(node, std::move(entries)).first;
            }

            // Strand: the first node at which a target appears decides it,
            // forward winning over reverse -- detect_strand's rule exactly.
            if (undecided > 0) {
                node_has_fwd.clear();
                for (const Entry& e : it->second)
                    for (uint32_t t : owners[e.pid])
                        if (strand[t] < 0) node_has_fwd[t] = node_has_fwd[t] || !e.rev;
                for (const auto& kv : node_has_fwd) {
                    strand[kv.first] = kv.second ? 0 : 1;
                    undecided--;
                }
            }

            // Forward-strand anchors. A target only has entries at or after the
            // node that decided it, so this never feeds an undecided target.
            for (const Entry& e : it->second) {
                if (e.rev) continue;
                for (uint32_t t : owners[e.pid]) {
                    if (strand[t] != 0) continue;
                    Anchor a;
                    a.read_begin = mappings[i].read_begin_offset;
                    a.read_end = mappings[i].read_end_offset;
                    a.path_offset = e.path_offset;
                    a.edge_node = node;
                    a.edge_offset = e.k;
                    a.source_mapping = i;
                    by_path[t][e.pid].push_back(a);
                }
            }
        }

        // ── flipped pass, shared by every reverse-strand target ──────────
        std::unordered_set<size_t> want_rev;
        for (size_t t = 0; t < T; t++)
            if (strand[t] == 1) want_rev.insert(targets[t].begin(), targets[t].end());
        if (!want_rev.empty()) {
            std::unordered_map<size_t, std::vector<Anchor>> rev_by_path;
            scan(reverse_strand_mappings(mappings), want_rev, rev_by_path, st);
            for (auto& kv : rev_by_path)
                for (uint32_t t : owners[kv.first])
                    if (strand[t] == 1) by_path[t][kv.first] = kv.second;
        }

        for (size_t t = 0; t < T; t++) {
            if (strand[t] < 0) continue;
            res[t].reverse_strand = (strand[t] == 1);
            res[t].touched_paths = by_path[t].size();
            finalize(by_path[t], res[t]);
            st.anchors += res[t].anchors.size();
        }
        st.ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        return res;
    }

private:
    /// Pick the path the read overlaps most and sort its anchors. Scores by
    /// DISTINCT covered mappings, not raw anchors, so a repeat-heavy path
    /// cannot out-score the one the read actually aligns to; ties go to the
    /// lowest path id so repeated runs agree.
    static void finalize(const std::unordered_map<size_t, std::vector<Anchor>>& by_path,
                         TargetResult& r) {
        const std::vector<Anchor>* best = nullptr;
        size_t best_pid = 0, best_cov = 0;
        for (const auto& kv : by_path) {
            std::unordered_set<size_t> covered;
            for (const Anchor& a : kv.second) covered.insert(a.source_mapping);
            if (best == nullptr || covered.size() > best_cov ||
                (covered.size() == best_cov && kv.first < best_pid)) {
                best = &kv.second; best_cov = covered.size(); best_pid = kv.first;
            }
        }
        if (best == nullptr) return;
        r.best_pid = best_pid;
        r.anchors = *best;
        std::sort(r.anchors.begin(), r.anchors.end(), [](const Anchor& a, const Anchor& b) {
            if (a.read_begin != b.read_begin) return a.read_begin < b.read_begin;
            return a.path_offset < b.path_offset;
        });
    }

    /// Which strand of the read aligns to a wanted path: 0 own orientation,
    /// 1 flipped, -1 neither. Returns on the first mapping that hits.
    /// One query per node settles both questions. Querying the node in the
    /// READ's orientation returns sequence 2p when the path traverses it the
    /// same way, and 2p+1 when it traverses it the other way -- so the
    /// orientation of the matching sequence IS the strand, and the flipped node
    /// never has to be queried.
    int detect_strand(const std::vector<SourceMapping>& mappings,
                      const std::unordered_set<size_t>& want, Stats& st) const {
        for (const SourceMapping& m : mappings) {
            const gbwt::node_type node = gbwt::Node::encode(
                static_cast<gbwt::size_type>(m.node_id), m.is_reverse);
            std::vector<gbwt::size_type> sa = sri_.decompressSA(node);
            st.sa_calls++;
            st.sa_entries += sa.size();
            int found = -1;
            for (gbwt::size_type v : sa) {
                const gbwt::size_type seq = sri_.seqId(v);
                if (!want.count(gbwt::Path::id(seq))) continue;
                if (!gbwt::Path::is_reverse(seq)) { found = 0; break; }  // forward wins
                found = 1;
            }
            if (found >= 0) return found;
        }
        return -1;
    }

    /// The one pass: decompressSA per distinct oriented node, cached because
    /// reads revisit nodes in segdups.
    void scan(const std::vector<SourceMapping>& mappings,
              const std::unordered_set<size_t>& want,
              std::unordered_map<size_t, std::vector<Anchor>>& by_path,
              Stats& st) const {
        struct Visit { size_t pid; size_t path_offset; size_t edge_offset; };
        std::unordered_map<gbwt::node_type, std::vector<Visit>> cache;

        for (size_t i = 0; i < mappings.size(); i++) {
            const gbwt::node_type node = gbwt::Node::encode(
                static_cast<gbwt::size_type>(mappings[i].node_id),
                mappings[i].is_reverse);

            auto it = cache.find(node);
            if (it == cache.end()) {
                std::vector<Visit> visits;
                std::vector<gbwt::size_type> sa = sri_.decompressSA(node);
                st.sa_calls++;
                st.sa_entries += sa.size();
                st.distinct_nodes++;
                for (size_t k = 0; k < sa.size(); k++) {
                    const gbwt::size_type seq = sri_.seqId(sa[k]);
                    if (gbwt::Path::is_reverse(seq)) continue;
                    const size_t pid = gbwt::Path::id(seq);
                    if (!want.count(pid)) continue;
                    // decompressSA measures from the sequence END; forwardStart
                    // turns that into a forward base position.
                    visits.push_back({pid,
                                      sri_.forwardStart(seq, sri_.seqOffset(sa[k]), node),
                                      k});
                }
                it = cache.emplace(node, std::move(visits)).first;
            }

            for (const Visit& v : it->second) {
                Anchor a;
                a.read_begin     = mappings[i].read_begin_offset;
                a.read_end       = mappings[i].read_end_offset;
                a.path_offset    = v.path_offset;
                a.edge_node      = node;
                a.edge_offset    = v.edge_offset;
                a.source_mapping = i;
                by_path[v.pid].push_back(a);
            }
        }
    }

    /// Identical to panindexer::make_reverse_strand_mappings.
    static std::vector<SourceMapping> reverse_strand_mappings(
        const std::vector<SourceMapping>& mappings) {
        size_t read_len = 0;
        for (const SourceMapping& m : mappings)
            read_len = std::max(read_len, m.read_end_offset);
        std::vector<SourceMapping> out;
        out.reserve(mappings.size());
        for (size_t k = mappings.size(); k-- > 0; ) {
            SourceMapping s = mappings[k];
            s.is_reverse = !s.is_reverse;
            const size_t rb = mappings[k].read_begin_offset;
            const size_t re = mappings[k].read_end_offset;
            s.read_begin_offset = (read_len >= re) ? (read_len - re) : 0;
            s.read_end_offset   = (read_len >= rb) ? (read_len - rb) : 0;
            out.push_back(s);
        }
        return out;
    }

    const gbwtgraph::GBZ& gbz_;
    const gbwt::SequenceLocate& sri_;
};

} // namespace srianchor
