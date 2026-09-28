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
        if (by_path.empty()) {
            st.ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t0).count();
            return {};
        }

        // Score by DISTINCT covered mappings, not raw anchors: a repeat-heavy
        // path would otherwise out-score the one the read actually aligns to.
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
        if (best_pid_out) *best_pid_out = best_pid;

        std::vector<Anchor> out = *best;
        std::sort(out.begin(), out.end(), [](const Anchor& a, const Anchor& b) {
            if (a.read_begin != b.read_begin) return a.read_begin < b.read_begin;
            return a.path_offset < b.path_offset;
        });
        st.anchors = out.size();
        st.ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t0).count();
        return out;
    }

private:
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
