#include <iostream>
#include <filesystem>
#include <fstream> // std::ofstream
#include <string>
#include <tuple>
#include <zlib.h>  
#include <chrono>
#include <cassert>

#include "command-line-parsing/founderblockgraph.h" // gengetopt-generated parser
#include "algo.hpp"
#include "utils.hpp"
#include "index.hpp"

#define FOUNDERBLOCKGRAPH_DEBUG

// GFAKluge setup
//#define GFAK_LINK_TYPE 1
//using gfak::GFAKluge;

using std::cerr, std::cout, std::flush, std::endl, std::flush, std::chrono::steady_clock, std::ofstream;
using std::vector, std::string, std::pair, std::tuple, std::tie, std::get, std::filesystem::path, std::ostringstream;
using fbg::utils::open_msa_file, fbg::utils::close_msa_file, fbg::utils::get_msa_line, fbg::utils::count_ones;
using fbg::algo::efg, fbg::algo::segmentation;
using fbg::index::msa_index;

int main(int argc, char* argv[])
{
	gengetopt_args_info argsinfo;
	if (cmdline_parser(argc, argv, &argsinfo) != 0) exit(1);

	if (argsinfo.non_elastic_flag) { cerr << "--non-elastic mode is not yet implemented!" << endl; exit(1); };
	if (argsinfo.gap_limit_arg != -1) { cerr << "gap limit is not yet implemented!" << endl; exit(1); };
	const string ignorechars = ((argsinfo.ignore_chars_arg == NULL) ? "" : string(argsinfo.ignore_chars_arg));
	if (argsinfo.heuristic_subset_arg == 0 or argsinfo.heuristic_subset_arg < -1) { cerr << "ERROR: wrong value of heuristic subset size!" << endl; exit(1); };
	if (argsinfo.threads_arg == 0 or argsinfo.threads_arg < -1) { cerr << "ERROR: wrong number of threads!" << endl; exit(1); };
	const path tmpdir(argsinfo.tmp_dir_arg);
	path in_file(argsinfo.input_arg);
	path out_file(argsinfo.output_arg);

	std::chrono::time_point<std::chrono::steady_clock> start, end;
	unsigned long long m, n;

	if (argsinfo.heuristic_subset_arg == -1) {
		// default elastic mode, single or multi-threaded
		cerr << "INFO: indexing the MSA..." << flush;
		start = steady_clock::now();
		msa_index index = fbg::index::index_external_memory(argsinfo.input_arg, tmpdir, ignorechars, m, n);
		end = steady_clock::now();
		cerr << " done. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;
		cerr << "INFO: indexed MSA[1.." << m << ", 1.." << n << "]" << endl; // TODO space usage

		cerr << "INFO: computing the optimal segmentation..." << flush;
		start = steady_clock::now();
		segmentation S = fbg::algo::optimal_segmentation_minmaxlength(m, n, index, argsinfo.threads_arg, (ignorechars != ""), argsinfo.disable_tricks_flag);
		end = steady_clock::now();
		cerr << " done. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;

		cerr << "INFO: computing and outputting graph..." << flush;
		start = steady_clock::now();
		efg g = fbg::algo::segment_msa(argsinfo.input_arg, n, S);
		ofstream out_stream(out_file, std::ios::out | std::ios::trunc);
		fbg::algo::output_msa_info(m, n, out_stream);
		fbg::algo::output_segmentation(S, out_stream);
		fbg::algo::output_block_info(g, out_stream);
		fbg::algo::output_efg(g, out_stream);
		if (argsinfo.output_paths_flag) fbg::algo::output_paths(in_file, S, g, out_stream);
		end = steady_clock::now();
		cerr << " done. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;
	} else {
		// heuristic subset mode, single or multi-threaded
		cerr << "INFO: computing a heuristic segmentation..." << flush;
		start = steady_clock::now();
		segmentation S = fbg::algo::heuristic_subset_segmentation_minmaxlength(in_file, tmpdir, argsinfo.heuristic_subset_arg, m, n, argsinfo.threads_arg, ignorechars, argsinfo.disable_tricks_flag);
		end = steady_clock::now();
		cerr << " done. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;

		cerr << "INFO: computing graph..." << flush;
		start = steady_clock::now();
		efg g = fbg::algo::segment_msa(argsinfo.input_arg, n, S);
		end = steady_clock::now();
		cerr << " done. (" << fbg::index::count_nodes(g) << " nodes, " << fbg::index::count_edges(g) << " edges, " << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;

		int iters = 0;
		cerr << "INFO: fixing the segmentation until indexable..." << flush;
		start = steady_clock::now();
		{
			vector<bool> to_remove;
			while (!fbg::algo::is_indexable(S, g, ignorechars, argsinfo.threads_arg, tmpdir, to_remove)) {
				cerr << " " << count_ones(to_remove) << " blocks to remove...";
				fbg::algo::prune(S, to_remove);
				g = fbg::algo::segment_msa(argsinfo.input_arg, n, S);
				cerr << fbg::index::count_nodes(g) << " nodes, " << fbg::index::count_edges(g) << " edges...";
				iters += 1;
			}
		}
		end = steady_clock::now();
		cerr << " done in " << iters << " iterations. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;

		cerr << "INFO: recomputing and outputting graph..." << flush;
		start = steady_clock::now();
		g = fbg::algo::segment_msa(argsinfo.input_arg, n, S);
		ofstream out_stream(out_file, std::ios::out | std::ios::trunc);
		fbg::algo::output_msa_info(m, n, out_stream);
		fbg::algo::output_segmentation(S, out_stream);
		fbg::algo::output_block_info(g, out_stream);
		fbg::algo::output_efg(g, out_stream);
		if (argsinfo.output_paths_flag) fbg::algo::output_paths(in_file, S, g, out_stream);
		end = steady_clock::now();
		cerr << " done. (" << std::chrono::duration<double>({ end - start }) << " seconds)." << endl;
	}

	return 0;
}
