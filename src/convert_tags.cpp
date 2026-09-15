#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

#include "pangenome_index/tag_arrays.hpp"
#include <gbwtgraph/gbz.h>

using namespace std;
using namespace panindexer;
using handlegraph::pos_t;

static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <input_algorithm_file> <output_compressed_file>"
         << " (--gbz <graph.gbz> | --num-seq N)\n"
         << "\n"
         << "  --gbz <graph.gbz>  derive the sequence count from the GBZ (preferred)\n"
         << "  --num-seq N        number of BWT sequences, i.e. endmarkers; for a\n"
         << "                     bidirectional GBWT this is 2 x the number of paths\n"
         << "\n"
         << "One of the two is REQUIRED. The BWT reserves one position per sequence\n"
         << "for its endmarker, and this tool prepends them. Omitting the count builds\n"
         << "a tag array (N-1) positions short of the r-index, which silently shifts\n"
         << "every tag and makes coordinate translation lose almost all of its bases.\n";
}

int main(int argc, char** argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const string input_file = argv[1];
    const string output_file = argv[2];

    // Default sidecar tmp paths; we always auto-generate
    string encoded_starts_tmp = output_file + ".encoded_starts.tmp";
    string bwt_intervals_tmp = output_file + ".bwt_intervals.tmp";

    // Flags parsing (after 2 required positional args)
    // Number of BWT sequences (= endmarkers) to prepend. There is no sane
    // default: 0 produces a tag array that is (sequences - 1) positions shorter
    // than the r-index it will be queried against, so every lookup past the
    // first endmarker reads a position belonging to some other sequence. That
    // failed silently -- as "gap" tags and dropped bases -- rather than loudly,
    // so the count is now required.
    size_t bwt_offset = 0;
    bool have_offset = false;
    string gbz_file_path;
    for (int i = 3; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "--gbz") {
            if (i + 1 >= argc) {
                cerr << "Error: missing value for " << arg << "\n";
                return 1;
            }
            gbz_file_path = argv[++i];
        } else if (arg == "--num-seq" || arg == "-n") {
            if (i + 1 >= argc) {
                cerr << "Error: missing value for " << arg << "\n";
                return 1;
            }
            try {
                bwt_offset = static_cast<size_t>(stoull(argv[++i]));
                have_offset = true;
            } catch (const std::exception&) {
                cerr << "Error: invalid value for --num-seq (expected unsigned integer)\n";
                return 1;
            }
        } else {
            cerr << "Error: unknown option '" << arg << "'\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (!gbz_file_path.empty()) {
        gbwtgraph::GBZ gbz;
        std::ifstream gbz_in(gbz_file_path, std::ios::binary);
        if (!gbz_in) {
            cerr << "Error: cannot open GBZ file: " << gbz_file_path << "\n";
            return 1;
        }
        try {
            gbz.simple_sds_load(gbz_in);
        } catch (const std::exception& e) {
            cerr << "Error: failed to load GBZ: " << e.what() << "\n";
            return 1;
        }
        const size_t from_gbz = gbz.index.sequences();
        if (have_offset && bwt_offset != from_gbz) {
            cerr << "Error: --num-seq " << bwt_offset << " disagrees with the GBZ, "
                 << "which has " << from_gbz << " sequences\n";
            return 1;
        }
        bwt_offset = from_gbz;
        have_offset = true;
        cerr << "Sequence count from GBZ: " << bwt_offset
             << " (endmarker positions to prepend)\n";
    }

    if (!have_offset) {
        cerr << "Error: the sequence count is required -- pass --gbz <graph.gbz> "
                "or --num-seq N.\n"
                "Without it the tag array is built (sequences - 1) positions "
                "shorter than the r-index,\nwhich shifts every tag and makes "
                "coordinate translation lose nearly all of its bases.\n";
        usage(argv[0]);
        return 1;
    }

    // Read the entire input file (algorithm.hpp format: raw gbwt::ByteCode runs) into memory
    vector<gbwt::byte_type> encoded_bytes;
    try {
        ifstream in(input_file, ios::binary);
        if (!in) {
            cerr << "Error: cannot open input file: " << input_file << "\n";
            return 1;
        }
        in.seekg(0, ios::end);
        streampos sz = in.tellg();
        if (sz < 0) {
            cerr << "Error: failed to stat input file size\n";
            return 1;
        }
        encoded_bytes.resize(static_cast<size_t>(sz));
        in.seekg(0, ios::beg);
        if (!encoded_bytes.empty()) {
            in.read(reinterpret_cast<char*>(encoded_bytes.data()), encoded_bytes.size());
            if (!in) {
                cerr << "Error: failed to read input file: " << input_file << "\n";
                return 1;
            }
        }
        in.close();
    } catch (const std::exception& e) {
        cerr << "Exception while reading input: " << e.what() << "\n";
        return 1;
    }

    // Open sidecar outputs (encoded_starts and bwt_intervals). Encoded runs go to a temp int_vector file.
    ofstream encoded_starts_out(encoded_starts_tmp, ios::binary | ios::trunc);
    if (!encoded_starts_out) {
        cerr << "Error: cannot open tmp file: " << encoded_starts_tmp << "\n";
        return 1;
    }

    ofstream bwt_intervals_out(bwt_intervals_tmp, ios::binary | ios::trunc);
    if (!bwt_intervals_out) {
        cerr << "Error: cannot open tmp file: " << bwt_intervals_tmp << "\n";
        return 1;
    }

    TagArray tag_array;
    
    // Pass 1: find max node id to determine width bits for int_vector encoding
    // Skip runs with node_id=0 (endmarkers) since we handle them separately
    uint64_t scan_loc = 0;
    uint64_t max_node_id = 0;
    while (scan_loc < encoded_bytes.size()) {
        gbwt::size_type decc = gbwt::ByteCode::read(encoded_bytes, scan_loc);
        auto decoded = TagArray::decode_run(decc);
        uint64_t nid = gbwtgraph::id(decoded.first);
        // Skip node_id=0 runs (endmarkers) when computing max_node_id
        // We handle endmarkers separately via --num-seq, so they shouldn't affect bit width
        if (nid > 0 && nid > max_node_id) max_node_id = nid;
    }
    // If max_node_id is still 0, we need at least 1 bit to represent node_id=0 for endmarkers
    if (max_node_id == 0 && bwt_offset > 0) {
        max_node_id = 0; // Keep at 0, compute_bits(0) returns 1 which is correct
    }
    auto compute_bits = [](uint64_t x) -> size_t {
        if (x == 0) return 1;
        return 64 - __builtin_clzll(x);
    };
    size_t node_bits = compute_bits(max_node_id);
    size_t width_bits = 10 + 1 + node_bits; // offset(10) + is_rev(1) + node id bits

    // Prepare temp file for encoded runs (sdsl int_vector<0>) and stream runs
    const std::string encoded_runs_path = output_file + ".encoded_runs.tmp";
    tag_array.begin_encoded_runs_sdsl(encoded_runs_path, width_bits);

    // Pass 2: stream all runs into encoded_runs_iv and sidecars
    // Note: runs are already merged by build_tags, BUT we need to check if the first run
    // should be merged with endmarkers if bwt_offset > 0
    
    // Track the previous run to handle merging with endmarkers
    pos_t previous_tag;
    uint16_t previous_length = 0;
    bool has_previous = false;
    
    // Prepend endmarkers if requested
    // Handle overflow: if bwt_offset > 65535, emit multiple runs of 65535 each
    if (bwt_offset > 0) {
        previous_tag = pos_t{0, 0, 0};
        size_t remaining = bwt_offset;
        while (remaining > 65535) {
            tag_array.append_compact_run_streamed(previous_tag, 65535, encoded_starts_out, bwt_intervals_out);
            remaining -= 65535;
        }
        previous_length = static_cast<uint16_t>(remaining);
        has_previous = true;
    }

    uint64_t loc = 0;
    size_t run_count = 0;
    while (loc < encoded_bytes.size()) {
        gbwt::size_type decc = gbwt::ByteCode::read(encoded_bytes, loc);
        auto decoded = TagArray::decode_run(decc);
        
        // Skip only pure endmarker runs (node_id=0, offset=0, is_rev=0); they are handled via --num-seq.
        // Runs with node_id=0 but non-zero offset/is_rev are kept as normal runs (e.g. gaps or special encoding).
        uint64_t nid = gbwtgraph::id(decoded.first);
        if (nid == 0) {
            std::cerr << "Skipping pure endmarker run: " << decoded.first << " " << decoded.second << std::endl;
            continue;  // pure endmarker, skip (we prepend endmarkers via --num-seq)
        }
        
        run_count++;
        
        // If current run has the same tag as previous, merge them (like merge_tags does)
        if (has_previous && previous_tag == decoded.first) {
            // Merge runs: add the lengths together
            uint32_t merged_length = static_cast<uint32_t>(previous_length) + static_cast<uint32_t>(decoded.second);
            // Handle uint16_t overflow
            while (merged_length > 65535) {
                tag_array.append_compact_run_streamed(previous_tag, 65535, encoded_starts_out, bwt_intervals_out);
                merged_length -= 65535;
            }
            previous_length = static_cast<uint16_t>(merged_length);
            // previous_tag stays the same
        } else {
            // Write previous run if it exists, then start new one
            if (has_previous) {
                tag_array.append_compact_run_streamed(previous_tag, previous_length, encoded_starts_out, bwt_intervals_out);
            }
            previous_tag = decoded.first;
            previous_length = decoded.second;
            has_previous = true;
        }
    }

    // Write the last run
    if (has_previous) {
        tag_array.append_compact_run_streamed(previous_tag, previous_length, encoded_starts_out, bwt_intervals_out);
    }

    // Finish encoded runs iv
    tag_array.end_encoded_runs_sdsl();
    encoded_starts_out.flush();
    encoded_starts_out.close();
    bwt_intervals_out.flush();
    bwt_intervals_out.close();

    // Write encoded runs iv into main index file, then append sdsl sidecars
    try {
        // Write encoded_runs_iv first
        {
            std::ofstream main_out(output_file, std::ios::binary | std::ios::trunc);
            if (!main_out) { cerr << "Error: cannot open output file: " << output_file << "\n"; return 1; }
            std::ifstream iv_in(encoded_runs_path, std::ios::binary);
            if (!iv_in) { cerr << "Error: cannot open temp encoded runs: " << encoded_runs_path << "\n"; return 1; }
            tag_array.load_encoded_runs_sdsl(iv_in);
            tag_array.serialize_encoded_runs_sdsl(main_out);
        }
        std::remove(encoded_runs_path.c_str());

        // Append encoded_starts_sd and bwt_intervals
        tag_array.merge_compressed_files_sdsl(output_file, encoded_starts_tmp, bwt_intervals_tmp);
    } catch (const std::exception& e) {
        cerr << "Error while merging compressed files: " << e.what() << "\n";
        return 1;
    }

    // Verification: Try to load the converted file
    cerr << "\n=== Verification: Testing file loading ===" << endl;
    try {
        TagArray verify_array;
        std::ifstream verify_in(output_file, std::ios::binary);
        if (!verify_in) {
            cerr << "ERROR: Cannot open converted file for verification: " << output_file << endl;
            return 1;
        }
        verify_array.load_compressed_tags_compact(verify_in);
        cerr << "✓ Successfully loaded converted file!" << endl;
        cerr << "  - Number of runs: " << verify_array.number_of_runs_compressed() << endl;
        cerr << "  - Bytes in encoded runs: " << verify_array.bytes_encoded_runs() << endl;
        cerr << "  - Bytes in encoded starts: " << verify_array.bytes_encoded_runs_starts_sd() << endl;
        cerr << "  - Bytes in BWT intervals: " << verify_array.bytes_bwt_intervals() << endl;
        cerr << "  - Total bytes: " << verify_array.bytes_total_compressed() << endl;
        verify_in.close();
    } catch (const std::exception& e) {
        cerr << "ERROR: Failed to load converted file: " << e.what() << endl;
        return 1;
    }

    // Comparison: Try to load and compare with merge_tags output
    const string merge_tags_file = "whole_genome_tag_array_compressed.tags";
    cerr << "\n=== Comparison: Checking against merge_tags output ===" << endl;
    try {
        std::ifstream merge_file(merge_tags_file, std::ios::binary);
        if (!merge_file) {
            cerr << "Note: Cannot find merge_tags file (" << merge_tags_file << ") for comparison" << endl;
        } else {
            TagArray merge_array;
            merge_array.load_compressed_tags_compact(merge_file);
            cerr << "✓ Successfully loaded merge_tags file!" << endl;
            cerr << "  - Number of runs: " << merge_array.number_of_runs_compressed() << endl;
            cerr << "  - Bytes in encoded runs: " << merge_array.bytes_encoded_runs() << endl;
            cerr << "  - Bytes in encoded starts: " << merge_array.bytes_encoded_runs_starts_sd() << endl;
            cerr << "  - Bytes in BWT intervals: " << merge_array.bytes_bwt_intervals() << endl;
            cerr << "  - Total bytes: " << merge_array.bytes_total_compressed() << endl;

            // Now compare with converted file
            TagArray converted_array;
            std::ifstream converted_in(output_file, std::ios::binary);
            converted_array.load_compressed_tags_compact(converted_in);
            converted_in.close();

            bool matches = true;
            if (converted_array.number_of_runs_compressed() != merge_array.number_of_runs_compressed()) {
                cerr << "✗ MISMATCH: Number of runs differs: " 
                     << converted_array.number_of_runs_compressed() << " vs " 
                     << merge_array.number_of_runs_compressed() << endl;
                matches = false;
            }
            if (converted_array.bytes_encoded_runs() != merge_array.bytes_encoded_runs()) {
                cerr << "✗ MISMATCH: Encoded runs bytes differ: " 
                     << converted_array.bytes_encoded_runs() << " vs " 
                     << merge_array.bytes_encoded_runs() << endl;
                matches = false;
            }
            if (converted_array.bytes_encoded_runs_starts_sd() != merge_array.bytes_encoded_runs_starts_sd()) {
                cerr << "✗ MISMATCH: Encoded starts bytes differ: " 
                     << converted_array.bytes_encoded_runs_starts_sd() << " vs " 
                     << merge_array.bytes_encoded_runs_starts_sd() << endl;
                matches = false;
            }
            if (converted_array.bytes_bwt_intervals() != merge_array.bytes_bwt_intervals()) {
                cerr << "✗ MISMATCH: BWT intervals bytes differ: " 
                     << converted_array.bytes_bwt_intervals() << " vs " 
                     << merge_array.bytes_bwt_intervals() << endl;
                matches = false;
            }

            // Compare actual run data using for_each_run_compact (always collect and compare)
            std::vector<std::pair<pos_t, uint64_t>> converted_runs, merge_runs;
            converted_array.for_each_run_compact([&](pos_t p, uint64_t len) {
                converted_runs.push_back({p, len});
            });
            merge_array.for_each_run_compact([&](pos_t p, uint64_t len) {
                merge_runs.push_back({p, len});
            });

            if (converted_runs.size() != merge_runs.size()) {
                cerr << "✗ MISMATCH: Run count from iteration differs: "
                     << converted_runs.size() << " vs " << merge_runs.size() << endl;
                matches = false;
            }

            // Compare runs at each index (up to min size) and print differing runs
            const size_t max_print = 100;
            size_t diff_count = 0;
            size_t compare_len = (converted_runs.size() < merge_runs.size())
                ? converted_runs.size() : merge_runs.size();
            for (size_t i = 0; i < compare_len; ++i) {
                const auto& conv = converted_runs[i];
                const auto& merg = merge_runs[i];
                if (conv.first != merg.first || conv.second != merg.second) {
                    matches = false;
                    diff_count++;
                    if (diff_count <= max_print) {
                        cerr << "✗ Run " << i << " differs:\n"
                             << "    converted: node_id=" << gbwtgraph::id(conv.first)
                             << " is_rev=" << gbwtgraph::is_rev(conv.first)
                             << " offset=" << gbwtgraph::offset(conv.first)
                             << " len=" << conv.second << "\n"
                             << "    merge:     node_id=" << gbwtgraph::id(merg.first)
                             << " is_rev=" << gbwtgraph::is_rev(merg.first)
                             << " offset=" << gbwtgraph::offset(merg.first)
                             << " len=" << merg.second << endl;
                    }
                }
            }
            if (diff_count > 0) {
                cerr << "  Total differing runs (in comparable range): " << diff_count;
                if (diff_count > max_print) {
                    cerr << " (printed first " << max_print << ")";
                }
                cerr << endl;
            }
            if (converted_runs.size() != merge_runs.size()) {
                size_t extra = (converted_runs.size() > merge_runs.size())
                    ? (converted_runs.size() - merge_runs.size()) : (merge_runs.size() - converted_runs.size());
                size_t end_idx = (converted_runs.size() > merge_runs.size() ? converted_runs.size() : merge_runs.size()) - 1;
                cerr << "  " << (converted_runs.size() > merge_runs.size() ? "Converted" : "Merge")
                     << " has " << extra << " extra run(s) (indices " << compare_len << ".." << end_idx << ")" << endl;
                // Print first few extra runs
                size_t from = compare_len;
                size_t to = (converted_runs.size() > merge_runs.size()) ? converted_runs.size() : merge_runs.size();
                for (size_t i = from; i < to && (i - from) < 5; ++i) {
                    const auto& r = (converted_runs.size() > merge_runs.size()) ? converted_runs[i] : merge_runs[i];
                    cerr << "    run[" << i << "]: node_id=" << gbwtgraph::id(r.first)
                         << " is_rev=" << gbwtgraph::is_rev(r.first)
                         << " offset=" << gbwtgraph::offset(r.first) << " len=" << r.second << endl;
                }
                if (to - from > 5) {
                    cerr << "    ... and " << (to - from - 5) << " more" << endl;
                }
            }

            if (matches) {
                cerr << "\n✓ Files match! Converted file is identical to merge_tags output." << endl;
            } else {
                cerr << "\n✗ Files differ! There are differences between the converted file and merge_tags output." << endl;
            }
            merge_file.close();
        }
    } catch (const std::exception& e) {
        cerr << "ERROR: Failed to load/compare with merge_tags file: " << e.what() << endl;
        return 1;
    }

    cerr << "\n=== Conversion complete ===" << endl;
    return 0;
}


