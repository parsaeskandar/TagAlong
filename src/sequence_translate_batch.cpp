// Batch runner for old-vs-new comparison.
//   sequence_translate_batch <graph.gbz> <graph.sri> <queries.tsv> <out.tsv> [--scan]
// queries.tsv: id \t src \t start \t end \t tgt
// Emits: id \t n_pairs \t xor64 \t sum64 \t ms \t contigs
//
// Per-base pairs are hashed ORDER-INDEPENDENTLY (xor and sum of a per-pair
// FNV-1a) so the two implementations can be compared without either having to
// emit or sort millions of lines, and without depending on output order.
#include "sequence_translation.cpp"
#include <fstream>
#include <sstream>
#include <set>

static inline uint64_t fnv1a(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return h;
}

int main(int argc, char** argv) {
    if (argc < 5) { std::cerr << "usage: " << argv[0]
        << " <graph.gbz> <graph.sri> <queries.tsv> <out.tsv> [--scan]\n"; return 1; }
    seqtrans::Mode mode = seqtrans::Mode::WALK;
    for (int i = 5; i < argc; i++) if (std::string(argv[i]) == "--scan") mode = seqtrans::Mode::SCAN;

    auto t = std::chrono::steady_clock::now();
    gbwtgraph::GBZ gbz; sdsl::simple_sds::load_from(gbz, argv[1]);
    gbwt::SequenceLocate sri; sdsl::load_from_file(sri, argv[2]);
    const gbwtgraph::GBWTGraph* g = &gbz.graph;
    sri.setGBWT(gbz.index, [g](gbwt::node_type n) -> gbwt::size_type {
        return g->get_length(gbwtgraph::GBWTGraph::node_to_handle(n)); });
    if (sri.sequence_length.size() != gbz.index.sequences() || sri.size() == 0) {
        std::cerr << "FATAL: .sri/GBZ mismatch\n"; return 2; }
    std::cerr << "[load] " << std::chrono::duration<double>(
        std::chrono::steady_clock::now()-t).count() << "s\n";

    seqtrans::Translator tr(gbz, sri);
    std::unordered_map<std::string, std::vector<gbwt::size_type>> tcache;
    std::ifstream in(argv[3]); std::ofstream out(argv[4]);
    out << "id\tn\txor64\tsum64\tms\tcontigs\n";
    std::string line; size_t done = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string id, src, tgt, a, b;
        std::getline(ss,id,'\t'); std::getline(ss,src,'\t');
        std::getline(ss,a,'\t');  std::getline(ss,b,'\t'); std::getline(ss,tgt,'\t');
        auto it = tcache.find(tgt);
        if (it == tcache.end()) it = tcache.emplace(tgt, tr.resolve(tgt)).first;

        seqtrans::Stats st;
        auto res = it->second.empty() ? std::vector<seqtrans::Block>()
                 : tr.translate(src, std::stoull(a), std::stoull(b), it->second, st, mode);
        auto folded = seqtrans::Translator::fold(res);

        uint64_t x = 0, sm = 0; size_t n = 0; std::set<std::string> ctgs;
        for (const auto& bl : folded) {
            ctgs.insert(bl.target_contig);
            for (size_t i = 0; i < bl.src_end - bl.src_start; i++) {
                size_t sp = bl.src_start + i;
                size_t tp = (bl.strand=='+') ? bl.tgt_start + i : bl.tgt_end - 1 - i;
                uint64_t h = fnv1a(bl.target_contig + ":" + std::to_string(sp)
                                   + ">" + std::to_string(tp) + bl.strand);
                x ^= h; sm += h; n++;
            }
        }
        std::string cs; for (const auto& c : ctgs) { if (!cs.empty()) cs += ","; cs += c; }
        out << id << "\t" << n << "\t" << x << "\t" << sm << "\t" << st.ms
            << "\t" << (cs.empty() ? "-" : cs) << "\n";
        if (++done % 50 == 0) { out.flush(); std::cerr << "  " << done << " done\n"; }
    }
    std::cerr << "[done] " << done << " queries\n";
    return 0;
}
