// Driver for the SequenceLocate-only coordinate translation.
//
//   sequence_translate <graph.gbz> <graph.sri> <src> <start> <end> <tgt> [--blocks] [--repeat N]
//
// src/tgt are "SAMPLE#HAP#CONTIG" or "SAMPLE#HAP"; start/end are 0-based
// half-open CONTIG coordinates.
#include "sequence_translation.cpp"

int main(int argc, char** argv) {
    if (argc < 7) {
        std::cerr << "usage: " << argv[0]
                  << " <graph.gbz> <graph.sri> <src> <start> <end> <tgt>"
                     " [--blocks] [--repeat N] [--quiet]\n";
        return 1;
    }
    const std::string gbz_path = argv[1], sri_path = argv[2];
    const std::string src = argv[3], tgt = argv[6];
    const size_t start = std::stoull(argv[4]), end = std::stoull(argv[5]);
    bool blocks = false, quiet = false, pairs = false; int repeat = 1;
    seqtrans::Mode mode = seqtrans::Mode::WALK;
    for (int i = 7; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--blocks") blocks = true;
        else if (a == "--quiet") quiet = true;
        else if (a == "--pairs") pairs = true;          // per-base, for diffing
        else if (a == "--scan") mode = seqtrans::Mode::SCAN;
        else if (a == "--verify") mode = seqtrans::Mode::VERIFY;
        else if (a == "--repeat" && i + 1 < argc) repeat = std::atoi(argv[++i]);
    }

    auto t = std::chrono::steady_clock::now();
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_path);
    double gbz_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t).count();

    t = std::chrono::steady_clock::now();
    gbwt::SequenceLocate sri;
    sdsl::load_from_file(sri, sri_path);
    const gbwtgraph::GBWTGraph* graph = &gbz.graph;
    gbwt::SequenceLocate::length_function len =
        [graph](gbwt::node_type node) -> gbwt::size_type {
            return graph->get_length(gbwtgraph::GBWTGraph::node_to_handle(node));
        };
    sri.setGBWT(gbz.index, len);     // REQUIRED after load
    double sri_ms = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - t).count();

    // The .sri carries no checksum, and a mismatched pair would silently give
    // plausible-looking wrong coordinates, so check the one invariant we have.
    if (sri.sequence_length.size() != gbz.index.sequences() || sri.size() == 0) {
        std::cerr << "FATAL: .sri does not match this GBZ ("
                  << sri.sequence_length.size() << " vs " << gbz.index.sequences()
                  << " sequences, size=" << sri.size() << ")\n";
        return 2;
    }
    std::cerr << "[load] gbz " << gbz_ms/1000 << "s + sri " << sri_ms/1000
              << "s = " << (gbz_ms+sri_ms)/1000 << "s   paths=" << gbz.index.metadata.paths()
              << " sequences=" << gbz.index.sequences() << "\n";

    seqtrans::Translator tr(gbz, sri);
    auto tgt_ids = tr.resolve(tgt);
    if (tgt_ids.empty()) { std::cerr << "no paths match target '" << tgt << "'\n"; return 3; }
    std::cerr << "[target] '" << tgt << "' -> " << tgt_ids.size() << " path(s)\n";

    seqtrans::Stats st;
    std::vector<seqtrans::Block> res;
    for (int r = 0; r < repeat; r++) {
        seqtrans::Stats s;
        res = tr.translate(src, start, end, tgt_ids, s, mode);
        if (r == repeat - 1) st = s;
        std::cerr << "[run " << (r+1) << "] " << s.ms << " ms\n";
    }

    auto folded = seqtrans::Translator::fold(res);
    std::cerr << "[stats] fragments=" << st.fragments
              << " nodes=" << st.nodes_walked
              << " locateSequence=" << st.locate_sequence
              << " locateForward=" << st.locate_forward
              << " visits=" << st.visits_enumerated
              << " blocks=" << st.blocks
              << " folded=" << folded.size()
              << " self_check_fail=" << st.self_check_fail
              << " twalks=" << st.target_walks
              << " tsteps=" << st.target_walk_steps
              << " verify=" << st.verify_mismatch << "/" << st.verify_checked
              << " ms=" << st.ms << "\n";

    if (pairs) {
        // One line per source base, so the output can be diffed against the
        // existing implementation's per-base correspondences.
        for (const auto& b : seqtrans::Translator::fold(res)) {
            for (size_t i = 0; i < b.src_end - b.src_start; i++) {
                size_t sp = b.src_start + i;
                size_t tp = (b.strand == '+') ? b.tgt_start + i
                                              : b.tgt_end - 1 - i;
                std::cout << sp << "\t" << tp << "\t" << b.strand << "\n";
            }
        }
        return 0;
    }
    if (!quiet) {
        const auto& show = blocks ? folded : res;
        size_t n = 0;
        for (const auto& b : show) {
            std::cout << b.target_contig << "\t" << b.src_start << "\t" << b.src_end
                      << "\t" << b.tgt_start << "\t" << b.tgt_end << "\t" << b.strand << "\n";
            if (++n >= 40 && !blocks) { std::cout << "... (" << show.size() << " total)\n"; break; }
        }
    }
    return 0;
}
