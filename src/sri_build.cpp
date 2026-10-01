// Build a .sri for a GBZ.  sri_build <graph.gbz> <out.sri>
#include <gbwt/sequence_locate.h>
#include <gbwtgraph/gbz.h>
#include <sdsl/simple_sds.hpp>
#include <iostream>
int main(int argc, char** argv) {
    if (argc < 3) { std::cerr << "usage: " << argv[0] << " <graph.gbz> <out.sri>\n"; return 1; }
    gbwtgraph::GBZ gbz;
    sdsl::simple_sds::load_from(gbz, argv[1]);
    const gbwtgraph::GBWTGraph* g = &gbz.graph;
    gbwt::SequenceLocate::length_function len =
        [g](gbwt::node_type n) -> gbwt::size_type {
            return g->get_length(gbwtgraph::GBWTGraph::node_to_handle(n)); };
    gbwt::SequenceLocate sri(gbz.index, len);
    sdsl::store_to_file(sri, argv[2]);
    std::cerr << "built " << argv[2] << ": " << gbz.index.sequences()
              << " sequences, size=" << sri.size() << "\n";
    return 0;
}
