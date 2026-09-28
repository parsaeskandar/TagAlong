// Driver for the SequenceLocate-only surjection anchor builder.
//
//   sri_anchors <graph.gbz> <graph.sri> <mappings.tsv> <anchors.tsv> [--repeat N]
//
// The mappings file comes from the same gaf_to_source_mappings() call that
// feeds the existing builder, so the A/B compares anchor searches, not GAF
// parsers. Format: "# case <id> <target>" then node_id/is_reverse/read_begin/
// read_end rows.
#include "sri_anchor_builder.cpp"

#include <sdsl/io.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "usage: " << argv[0]
                  << " <graph.gbz> <graph.sri> <mappings.tsv> <anchors.tsv>"
                     " [--repeat N]\n";
        return 1;
    }
    const std::string gbz_path = argv[1], sri_path = argv[2];
    const std::string map_path = argv[3], out_path = argv[4];
    int repeat = 1;
    for (int i = 5; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--repeat" && i + 1 < argc) repeat = std::atoi(argv[++i]);
    }

    auto t = std::chrono::steady_clock::now();
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, gbz_path);
    const double gbz_s = std::chrono::duration<double>(
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
    const double sri_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t).count();

    // The .sri carries no checksum, so check the one invariant we have.
    if (sri.sequence_length.size() != gbz.index.sequences() || sri.size() == 0) {
        std::cerr << "FATAL: .sri does not match this GBZ ("
                  << sri.sequence_length.size() << " vs " << gbz.index.sequences()
                  << " sequences, size=" << sri.size() << ")\n";
        return 2;
    }
    std::cerr << "[load] gbz " << gbz_s << "s + sri " << sri_s << "s = "
              << (gbz_s + sri_s) << "s\n";

    srianchor::Builder builder(gbz, sri);

    std::ifstream in(map_path);
    if (!in) { std::cerr << "cannot read " << map_path << "\n"; return 3; }
    std::ofstream out(out_path);
    if (!out) { std::cerr << "cannot write " << out_path << "\n"; return 3; }
    out << "case\tread_begin\tread_end\tpath_offset\tedge_node\tedge_offset\n";

    std::string case_id, target;
    std::vector<srianchor::SourceMapping> mappings;

    auto flush_case = [&]() {
        if (case_id.empty()) return;
        std::vector<size_t> pids = builder.resolve(target);
        if (pids.empty()) {
            std::cerr << case_id << "\t" << target << "\tunknown_path\n";
            mappings.clear();
            return;
        }
        srianchor::Stats st;
        std::vector<srianchor::Anchor> anchors;
        size_t best_pid = 0;
        double best_ms = 0.0;
        for (int r = 0; r < repeat; r++) {
            srianchor::Stats s;
            anchors = builder.build(mappings, pids, s, &best_pid);
            // Fastest run: repeats remove first-touch page cache effects.
            if (r == 0 || s.ms < best_ms) best_ms = s.ms;
            st = s;
        }
        for (const auto& a : anchors) {
            out << case_id << '\t' << a.read_begin << '\t' << a.read_end << '\t'
                << a.path_offset << '\t' << a.edge_node << '\t'
                << a.edge_offset << '\n';
        }
        std::cerr << case_id << '\t' << target
                  << "\tanchors=" << anchors.size()
                  << "\tms=" << best_ms
                  << "\tresolved=" << pids.size()
                  << "\ttouched=" << st.touched_paths
                  << "\tbest=" << (anchors.empty() ? "-" : builder.path_name(best_pid))
                  << "\tsa_calls=" << st.sa_calls
                  << "\tsa_entries=" << st.sa_entries
                  << "\tnodes=" << st.mappings
                  << "\tdistinct=" << st.distinct_nodes
                  << "\trev=" << (st.reverse_strand ? 1 : 0)
                  << '\n';
        mappings.clear();
    };

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            flush_case();
            std::istringstream ss(line);
            std::string hash, kw;
            ss >> hash >> kw >> case_id >> target;
            continue;
        }
        std::istringstream ss(line);
        srianchor::SourceMapping m;
        int rev = 0;
        ss >> m.node_id >> rev >> m.read_begin_offset >> m.read_end_offset;
        m.is_reverse = (rev != 0);
        mappings.push_back(m);
    }
    flush_case();
    out.close();
    std::cerr << "[wrote] " << out_path << "\n";
    return 0;
}
