// 10k-query benchmark for the SequenceLocate-only translator.
#include <iomanip>
//   sequence_translate_bench <graph.gbz> <graph.sri> <N> [out.tsv]
// 90% of sources are GRCh38#0 / CHM13#0 contigs (the browser's real pattern),
// 10% are other haplotypes. Mixed lengths. One load, N queries.
#include "sequence_translation.cpp"
#include <random>
#include <fstream>
#include <map>

int main(int argc, char** argv) {
    if (argc < 4) { std::cerr << "usage: " << argv[0]
        << " <graph.gbz> <graph.sri> <N> [out.tsv]\n"; return 1; }
    const size_t N = std::stoull(argv[3]);
    const std::string outpath = (argc > 4) ? argv[4] : "";

    auto t = std::chrono::steady_clock::now();
    gbwtgraph::GBZ gbz; sdsl::simple_sds::load_from(gbz, argv[1]);
    gbwt::SequenceLocate sri; sdsl::load_from_file(sri, argv[2]);
    const gbwtgraph::GBWTGraph* graph = &gbz.graph;
    sri.setGBWT(gbz.index, [graph](gbwt::node_type n) -> gbwt::size_type {
        return graph->get_length(gbwtgraph::GBWTGraph::node_to_handle(n)); });
    if (sri.sequence_length.size() != gbz.index.sequences() || sri.size() == 0) {
        std::cerr << "FATAL: .sri/GBZ mismatch\n"; return 2; }
    std::cerr << "[load] " << std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t).count() << "s\n";

    // Collect contig names and their extents, and the haplotype names.
    const gbwt::Metadata& md = gbz.index.metadata;
    std::map<std::string, size_t> contig_len;            // sample#hap#contig -> max end
    std::vector<std::string> haps_all;
    { std::unordered_set<std::string> seen;
      for (gbwt::size_type pid = 0; pid < md.paths(); pid++) {
          gbwt::FullPathName fp = md.fullPath(pid);
          std::string ctg = fp.sample_name + "#" + std::to_string(fp.haplotype)
                          + "#" + fp.contig_name;
          size_t e = fp.offset + sri.sequenceLength(gbwt::Path::encode(pid, false));
          auto it = contig_len.find(ctg);
          if (it == contig_len.end() || it->second < e) contig_len[ctg] = e;
          std::string hap = fp.sample_name + "#" + std::to_string(fp.haplotype);
          if (seen.insert(hap).second) haps_all.push_back(hap);
      } }
    std::vector<std::string> ref_ctgs, other_ctgs;
    for (const auto& kv : contig_len) {
        if (kv.second < 2000000) continue;               // skip tiny contigs
        if (kv.first.rfind("GRCh38#0#",0)==0 || kv.first.rfind("CHM13#0#",0)==0)
             ref_ctgs.push_back(kv.first);
        else other_ctgs.push_back(kv.first);
    }
    std::cerr << "[setup] " << ref_ctgs.size() << " reference contigs, "
              << other_ctgs.size() << " other, " << haps_all.size() << " haplotypes\n";

    seqtrans::Translator tr(gbz, sri);
    std::mt19937_64 rng(20260925);
    const size_t LENS[] = {1000, 10000, 100000, 1000000};
    const int    WTS[]  = {30, 30, 30, 10};

    std::ofstream out;
    if (!outpath.empty()) { out.open(outpath); out << "i\tsrc\ttgt\tlen\tms\tblocks\tfolded\tnodes\tselffail\n"; }

    std::map<size_t, std::vector<double>> by_len;
    std::vector<double> all; all.reserve(N);
    size_t selffail = 0, empty = 0;
    // Target path ids are resolved per haplotype and cached: resolving scans
    // all 94k paths, which would otherwise dominate the benchmark.
    std::unordered_map<std::string, std::vector<gbwt::size_type>> tcache;

    auto t0 = std::chrono::steady_clock::now();
    for (size_t i = 0; i < N; i++) {
        bool use_ref = (rng() % 100) < 90 || other_ctgs.empty();
        const std::string& src = use_ref ? ref_ctgs[rng() % ref_ctgs.size()]
                                         : other_ctgs[rng() % other_ctgs.size()];
        const std::string& tgt = haps_all[rng() % haps_all.size()];
        int w = rng() % 100, acc = 0; size_t L = LENS[0];
        for (int k = 0; k < 4; k++) { acc += WTS[k]; if (w < acc) { L = LENS[k]; break; } }
        size_t clen = contig_len[src];
        size_t start = (clen > L + 1) ? (rng() % (clen - L - 1)) : 0;

        auto it = tcache.find(tgt);
        if (it == tcache.end()) it = tcache.emplace(tgt, tr.resolve(tgt)).first;
        if (it->second.empty()) { i--; continue; }

        seqtrans::Stats st;
        auto res = tr.translate(src, start, start + L, it->second, st);
        auto folded = seqtrans::Translator::fold(res);
        all.push_back(st.ms); by_len[L].push_back(st.ms);
        selffail += st.self_check_fail;
        if (res.empty()) empty++;
        if (out.is_open())
            out << i << "\t" << src << "\t" << tgt << "\t" << L << "\t" << st.ms
                << "\t" << st.blocks << "\t" << folded.size() << "\t"
                << st.nodes_walked << "\t" << st.self_check_fail << "\n";
        if ((i+1) % 1000 == 0) {
            auto s = all; std::sort(s.begin(), s.end());
            std::cerr << "  [" << (i+1) << "/" << N << "] "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count()
                      << "s  median " << s[s.size()/2] << "ms  p99 "
                      << s[(size_t)(s.size()*0.99)] << "ms  max " << s.back() << "ms\n";
        }
    }
    double total = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();

    auto pct = [](std::vector<double> v, double p) {
        std::sort(v.begin(), v.end());
        return v[std::min(v.size()-1, (size_t)(v.size()*p))]; };
    std::cerr << "\n===== " << N << " queries in " << total << "s =====\n";
    std::cerr << "   length      n    median       p90       p99       max\n";
    for (auto& kv : by_len)
        std::cerr << std::setw(9) << kv.first << std::setw(7) << kv.second.size()
                  << std::setw(10) << pct(kv.second,0.5) << std::setw(10) << pct(kv.second,0.9)
                  << std::setw(10) << pct(kv.second,0.99) << std::setw(10) << pct(kv.second,1.0) << "\n";
    std::cerr << std::setw(9) << "ALL" << std::setw(7) << all.size()
              << std::setw(10) << pct(all,0.5) << std::setw(10) << pct(all,0.9)
              << std::setw(10) << pct(all,0.99) << std::setw(10) << pct(all,1.0) << "\n";
    std::cerr << "  self_check_fail: " << selffail << "   empty results: " << empty << "\n";
    return 0;
}
