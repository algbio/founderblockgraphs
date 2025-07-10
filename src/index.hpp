#ifndef INDEX_HPP
#define INDEX_HPP
//#define INDEX_HPP_DEBUG

#include <filesystem>
#include <iostream>
#include <vector> // std::erase
#include <fstream> // ofstream
#include <tuple> // pair
#include <unordered_map>
#include <unordered_set>
#include <mutex> // mutex, scoped_lock

#include <sdsl/construct.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/util.hpp>

#include "utils.hpp"

namespace fbg::index {

using std::string, std::cerr, std::endl, std::vector, std::filesystem::path, std::filesystem::is_directory, std::filesystem::remove, std::ostringstream, std::ofstream, std::pair, std::unordered_map, std::unordered_set, std::to_string, std::scoped_lock, std::mutex;
using fbg::utils::open_msa_file, fbg::utils::close_msa_file, fbg::utils::get_msa_line, fbg::utils::get_rows, fbg::utils::get_path;
using sdsl::bit_vector, sdsl::rank_support_v5, sdsl::select_support_mcl;
typedef sdsl::cst_sct3<> cst_type;
typedef sdsl::csa_wt<> sa_type;
typedef sdsl::rank_support_v5<> rank_type;
typedef sdsl::select_support_mcl<> select_type;
typedef cst_type::node_type node_t;
typedef cst_type::size_type size_type;
typedef std::vector<unsigned long long>::size_type size_type;

struct msa_index {
	cst_type cst; // compressed suffix tree of concat(enation)
	bit_vector concat_separators; // bit vector marking the separator characters of concat
	vector<bit_vector> nongaps; // bit matrix marking non-gaps (1s) and gaps (0s)
	vector<bit_vector> ignores; // bit matrix marking ignore characters (1s?)
	rank_type rs_concat_separators; // rank support
	vector<rank_type> rs_nongaps, rs_ignores; // rank support
	vector<select_type> ss_nongaps, ss_ignores; // select support
};

/* requires: msa path, tmp dir, ignore chars
 * modifies: optional m, n to store msa size
 * returns: full msa index
 * notes: thread safe if cap != -1, stop after reading cap rows if cap != -1, TODO store/load index from memory? */
unsigned long long _m, _n, _startrow; // dummy variables
mutex msa_batch_m; // lock for reading contiguous MSA rows
msa_index index_external_memory(const path &msapath, const path &tmpdir, const string ignorechars, unsigned long long &m = _m, unsigned long long &n = _n, unsigned long long &startrow = _startrow, const int cap = -1)
{
	assert(is_directory(tmpdir));
	assert((cap == -1 and msapath != "") or cap > 0);

	msa_index index; // result
	path concat_file;
	{
		scoped_lock l(msa_batch_m);
		startrow = get_rows();
		concat_file = path((tmpdir / get_path().stem()).string() + ((cap != -1) ? to_string(startrow) + "-" + to_string(startrow + cap - 1) : "") + ".plain"); // row concat
		ofstream concat_stream(concat_file, std::ios::out | std::ios::trunc);
		if (cap == -1)
			open_msa_file(msapath);
		//scoped_lock l(msa_batch_m);
		m = 0;
		string line, _;
		if (!get_msa_line(line, _)) {
			remove(concat_file);
			return index; // empty index, but m = 0
		}
		n = line.size();
		do {
			bit_vector nongap(n, 1), ignore(n, 0);
			for (size_type j = 0; j < line.size(); j++) {
				if (line[j] == GAP_CHARACTER)
					nongap[j] = 0;
				if (ignorechars.find(line[j]) != std::string::npos)
					ignore[j] = 1;
			}
			index.nongaps.push_back(std::move(nongap));
			index.ignores.push_back(std::move(ignore));

			std::erase(line, GAP_CHARACTER);
			concat_stream << line << SEPARATOR_CHARACTER;
			const size_type count_previous = index.concat_separators.size();
			index.concat_separators.resize(index.concat_separators.size() + line.size() + 1);
			for (size_type k = 0; k < line.size(); k++)
				index.concat_separators[count_previous + k] = 0;
			index.concat_separators[index.concat_separators.size() - 1] = 1;
			m += 1;
		} while ((cap == -1 or m < (unsigned long long)cap) and get_msa_line(line, _));
		if (cap == -1)
			close_msa_file();
		concat_stream.close();
	}

	// build cst and rank/support structures
	sdsl::construct(index.cst, concat_file.string(), 1); // generate index
	index.rs_concat_separators = rank_type(&index.concat_separators);
	index.rs_nongaps.resize(m);
	for (size_type i = 0; i < m; i++)
		index.rs_nongaps[i] = rank_type(&index.nongaps[i]);
	index.ss_nongaps.resize(m);
	for (size_type i = 0; i < m; i++)
		index.ss_nongaps[i] = select_type(&index.nongaps[i]);
	index.rs_ignores.resize(m);
	for (size_type i = 0; i < m; i++)
		index.rs_ignores[i] = rank_type(&index.ignores[i]);
	index.ss_ignores.resize(m);
	for (size_type i = 0; i < m; i++)
		index.ss_ignores[i] = select_type(&index.ignores[i]);

	remove(concat_file);
	return index;
}

/* elastic founder/block graph with nodes partitioned into blocks and arbitrary edges, actually (so a layered DAG) */
struct efg {
	vector<unordered_map<string,unsigned long>> blocks;
	unordered_map<unsigned long,size_type> node_to_block; // node -> its block
	unordered_map<unsigned long,unordered_set<unsigned long>> adjacency_lists; // node -> out-neighbors
};

unsigned long count_nodes(const efg &g)
{
	unsigned long n = 0;
	for (const auto &b : g.blocks) n += b.size();
	return n;
}

unsigned long count_edges(const efg &g)
{
	unsigned long e = 0;
	for (const auto &[_,l] : g.adjacency_lists) e += l.size();
	return e;
}

/* index for matching and validation of indexability */
struct efg_index {
	vector<string> labels;
	vector<unsigned long> edge_targets;
	vector<bool> is_source, is_sink;
	sa_type concat; // concatenation of edge labels
	bit_vector concat_leaders; // 1s mark the start of a node's neighborhood
	bit_vector concat_separators; // 1s mark the delimiters
	// rank/select support
	rank_type rs_concat_leaders, rs_concat_separators;
	select_type ss_concat_separators;
};

/* modifies: removes from g.blocks the node labels and stores them in .labels */
efg_index index_efg_external_memory(efg &g, const path &tmpdir) {
	assert(is_directory(tmpdir));

	efg_index index; // result
	path concat_file((tmpdir / "graph_edge_concat.plain").string());
	ofstream concat_stream(concat_file, std::ios::out | std::ios::trunc);

	const unsigned long n = count_nodes(g);
	index.labels = vector<string>(n);
	index.is_source = vector<bool>(n, true);
	index.is_sink   = vector<bool>(n, true);

	for (size_type b = 0; b < g.blocks.size(); b++) {
		for (const auto &[label,node] : g.blocks[b]) {
			index.labels[node] = std::move(label);
		}
	}

	// edge index is concatenation of edge labels with separator char
	//for (size_type b = 0; b < g.blocks.size(); b++) {
	//	for (const auto &[label,node] : g.blocks[b]) {
	for (unsigned long node = 0; node < n; node++) {
		for (const auto out_neighbor : g.adjacency_lists[node]) {
			index.is_source[node] = false;
			index.is_sink[out_neighbor] = false;
			index.edge_targets.push_back(out_neighbor);

			concat_stream << index.labels[node] << index.labels[out_neighbor] << SEPARATOR_CHARACTER;
			const size_type previous_size = index.concat_separators.size();
			const size_type added = index.labels[node].size() + index.labels[out_neighbor].size() + 1;
			index.concat_separators.resize(previous_size + added);
			for (size_type k = 0; k < added - 1; k++) index.concat_separators[previous_size + k] = 0;
			index.concat_separators[previous_size + added - 1] = 1;
		}
		const size_type previous_size = index.concat_leaders.size();
		const size_type added = index.concat_separators.size() - previous_size;
		if (added > 0) {
			index.concat_leaders.resize(previous_size + added);
			index.concat_leaders[previous_size] = 1;
			for (size_type k = 1; k < added; k++) index.concat_leaders[previous_size + k] = 0;
		} else {
			concat_stream << SEPARATOR_CHARACTER;
			index.concat_separators.resize(previous_size + 1);
			index.concat_separators[previous_size-1] = 0;
			index.concat_separators[previous_size] = 1;
			index.concat_leaders.resize(previous_size + 1);
			index.concat_leaders[previous_size] = 1;
		}
	}
	concat_stream.close();

	// build cst and rank/support structures
	sdsl::construct(index.concat, concat_file.string(), 1); // generate index
	index.rs_concat_leaders = rank_type(&index.concat_leaders);
	index.rs_concat_separators = rank_type(&index.concat_separators);
	index.ss_concat_separators = select_type(&index.concat_separators);

#ifdef INDEX_HPP_DEBUG
	cerr << "\nnode labels are"; for (auto &l : index.labels) cerr << " " << l; cerr << endl;
	cerr << "edge targets are"; for (auto &e : index.edge_targets) cerr << " " << e; cerr << endl;
	cerr << "concat is            "; for (size_type i = 0; i < index.concat.size(); i++) cerr << index.concat.bwt[index.concat.isa[(i+1)%index.concat.size()]]; cerr << endl;
	cerr << "concat_leaders is    " << index.concat_leaders << endl;
	cerr << "concat_separators is " << index.concat_separators << endl;
#endif

	remove(concat_file);
	return index;
}

} // namespace fbg::index
#endif // ifndef INDEX_HPP
