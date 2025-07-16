#ifndef ALGO_HPP
#define ALGO_HPP
//#define ALGO_HPP_DEBUG

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <tuple>
#include <limits>
#include <list>
#include <algorithm> // sort, reverse
#include <filesystem>
#include <thread>
#include <functional> // ref
#include <atomic> // atomic

#include "sdsl/int_vector.hpp"
#include "utils.hpp"
#include "index.hpp"

namespace fbg::algo {

using std::cerr, std::endl, std::string, std::vector, std::unordered_map, std::ofstream, std::unordered_set, std::pair, std::list, std::max, std::min, std::filesystem::path, std::filesystem::exists, std::filesystem::create_directory, std::filesystem::remove, std::thread, std::jthread, std::ref, std::atomic;
using fbg::utils::open_msa_file, fbg::utils::get_msa_line, fbg::utils::close_msa_file, fbg::utils::contains_chars;
using fbg::index::msa_index, fbg::index::efg, fbg::index::efg_index, fbg::index::node_t, fbg::index::sa_type, fbg::index::index_external_memory;
using sdsl::bit_vector;

typedef unsigned long long seg_index;
const unsigned long long SEG_INDEX_MAX = std::numeric_limits<seg_index>::max();
typedef vector<unsigned long long> segmentation; // sorted vector starting at 0 and ending at n
typedef std::vector<unsigned long long>::size_type seg_size_t;
typedef std::vector<unsigned long long>::size_type size_type;

/* requires: segmentation S is sorted vector starting at 0 and ending at n
 * returns: elastic block graph (or a layered DAG if a segments contains the empty string)
 * notes: MSA is streamed from disk, graph (no paths) is kept in memory */
efg segment_msa(const path &msa_path, const unsigned long long n, const segmentation &S)
{
	assert(S.at(0) == 0 and S.back() == n);
#ifdef ALGO_HPP_DEBUG
	cerr << "DEBUG: segmentation is";
	for (auto &s : S) cerr << " " << s;
	cerr << endl;
#endif

	unsigned long nodes = 0;
	vector<unordered_map<string,unsigned long>> blocks(S.size() - 1);
	unordered_map<unsigned long,size_type> node_to_block;
	unordered_map<unsigned long,unordered_set<unsigned long>> adjacency_lists;

	open_msa_file(msa_path);
	string line, _dummy;
	while (get_msa_line(line, _dummy)) {
		assert(line.size() == n);

		seg_index prev = SEG_INDEX_MAX;
		for (seg_size_t i = 0; i < S.size() - 1; i++) {
			assert(S[i] < S[i+1]);

			string label = line.substr(S[i], S[i+1] - S[i]);
			std::erase(label, GAP_CHARACTER); // remove gaps

			if (label == "")
				continue;

			if (blocks[i].contains(label)) {
				// node label in this block is already present
				const unsigned long id = blocks[i][label];
				if (prev != SEG_INDEX_MAX) {
					assert(id != 0);
					adjacency_lists[prev].insert(id);
				}
				prev = id;
			} else {
				// new node
				const unsigned long newid = nodes++;
				blocks[i].insert({ label, newid });
				node_to_block.insert({ newid, i });
				adjacency_lists.insert({ newid, unordered_set<unsigned long>() });
				if (prev != SEG_INDEX_MAX) {
					assert(newid != 0);
					adjacency_lists[prev].insert(newid);
				}
				prev = newid;
			}
		}
	}
	close_msa_file();

#ifdef ALGO_HPP_DEBUG
	cerr << "DEBUG: blocks are ";
	for (auto &b : blocks) {
		cerr << "{";
		for (auto &[label,node] : b) {
			cerr << " " << node << ":" << label;
		}
		cerr << " }";
	}
	cerr << endl;
#endif

	return efg({ std::move(blocks), std::move(node_to_block), std::move(adjacency_lists) });
}

void output_msa_info(const unsigned long long m, const unsigned long long n, ofstream &out)
{
	out << "M\t" << m << "\t" << n << "\n";
}

void output_segmentation(const segmentation &S, ofstream &out)
{
	// 0-indexed to 1-indexed, only starting cols (see xGFAspec.md)
	out << "X";
	for (seg_size_t i = 0; i < S.size() - 1; i++)
		out << "\t" << S[i]+1;
	out << "\n";
}

void output_block_info(const efg &g, ofstream &out)
{
	out << "B";
	for (auto &b : g.blocks)
		out << "\t" << b.size();
	out << "\n";
}

/* TODO: rename vertices? */
void output_efg(const efg &g, ofstream &out)
{
	for (auto &b : g.blocks) {
		for (auto &[label, node] : b) {
			out << "S\t" << node << "\t" << label << "\n";
			for (const auto &outneighbor : g.adjacency_lists.at(node)) {
				out << "L\t" << node << "\t+\t" << outneighbor << "\t+\t0M" << "\n";
			}
		}
	}
}

/* requires: MSA not to have empty lines, TODO check this somewhere
 * notes: fully streams whole MSA from disk */
void output_paths(path &msa_path, const segmentation &S, const efg &g, ofstream &out)
{
	open_msa_file(msa_path);
	string line, id;
	while (get_msa_line(line, id)) {
		out << "P\t" << id << "\t";
		bool first = true;
		for (seg_size_t i = 0; i < S.size() - 1; i++) {
			string label = line.substr(S[i], S[i+1] - S[i]);
			std::erase(label, GAP_CHARACTER); // remove gaps
			if (label == "")
				continue;

			//assert(g.blocks[i].contains(label));
			const unsigned long node = g.blocks[i].at(label);

			out << ((first) ? "" : ",") << node << "+";
			first = false;
		}
		out << "\t*\n";
	}
}

void update_max(atomic<size_type> &x, const size_type value)
{
	// see stackoverflow.com/questions/16190078/how-to-atomically-update-a-maximum-value
	size_type prev_value = x;
	while (prev_value < value and !x.compare_exchange_weak(prev_value, value)) {};
}
void update_max(size_type &x, const size_type value)
{
	x = max(x, value);
}

/* requires: index to be generated with function index_external_memory
 * notes: T alternative is atomic<size_type> for multi-thread setting
 *        based on Algorithm 4 from 10.1016/j.tcs.2023.114269 plus tricks from 10.1093/bioinformatics/btaf225 */
template<typename T = size_type>
void compute_f(
		size_type const m,
		size_type const n,
		const msa_index &index,
		vector<T> &f, // store result here
		const bool ignorechars = false,
		const bool disable_efg_tricks = false
) {
	// unpacking index
	const auto &cst = index.cst;
	const auto &rs_concat_separators = index.rs_concat_separators;
	const auto &nongaps = index.nongaps;
	const auto &rs_nongaps = index.rs_nongaps;
	const auto &ss_nongaps = index.ss_nongaps;
	const auto &ignores = index.ignores;
	const auto &rs_ignores = index.rs_ignores;
	const auto &ss_ignores = index.ss_ignores;

	// find the nodes corresponding to reading each whole row
	vector<node_t> leaves(m, cst.root());
	unordered_map<size_type, size_type> leavesmap; // leaf index -> MSA row
	for (size_type next = 0, i = 0; i < m; i++) {
		leaves[i] = cst.select_leaf(cst.csa.isa[next] + 1);
		leavesmap[cst.lb(leaves[i])] = i;
		next += rs_nongaps[i].rank(n) + 1;
	}

	// binary coloring of the leaves
	bit_vector color(cst.size(cst.root()), false);
	bit_vector fullrow;
	if (disable_efg_tricks)
		fullrow = bit_vector(m, false); // mark if leaves[i] is still the initial value
	else
		fullrow = bit_vector(m, true); // mark if leaves[i] is still the initial value

	for (size_type x = 0; x < n; x++) {
#ifdef ALGO_HPP_DEBUG
		cerr << "DEBUG: computing f[x], with x = " << x << endl;
		cerr << "DEBUG: leaves for f[" << x << "] are: ";
		for (auto l : leaves) cerr << l << " ";
		cerr << endl;
#endif

		size_type fimax = x;
		// Mark each leaf in leaves, after filtering
		for (size_type i = 0; i < m; i++) {
			if (!disable_efg_tricks and fullrow[i]) // leaf corresponds to full row
				continue;
			if (ignorechars and ignores[i][x]) // leaf corresponds to suffix starting with ignorechar
				continue;
			for (size_type ll = cst.lb(leaves[i]); ll <= cst.rb(leaves[i]); ll++) { // cst is not a generalized suffix tree
				color[ll] = true;
			}
		}

		// Process each set of contiguous leaves
		for (size_type i = 0; i < m; i++) {
			node_t const l = leaves[i];
			if (fullrow[i] or (ignorechars and ignores[i][x])) // leaves corresponding to full rows or ignorechar
				continue;
			if (cst.lb(l) == 0 || color[cst.lb(l) - 1] == false) {
				// if leftmost leaf does not correspond to row i, skip
				if (rs_concat_separators.rank(cst.sn(cst.select_leaf(cst.lb(l) + 1))) != i)
					break;
				size_type lb = cst.lb(l);
				size_type rb = cst.rb(l);
				while (rb < cst.size(cst.root()) - 1 && color[cst.rb(cst.select_leaf(rb)) + 2]) {
					rb = rb + 1;
				}

				// Find the exclusive ancestors
				node_t w = l;
				while (cst.rb(w) <= rb) {
					node_t parent = cst.parent(w);
					if (lb <= cst.lb(parent) && cst.rb(parent) <= rb) {
						// parent is a correct replacement
						w = parent;
					} else {
						// parent fails so w is an exclusive ancestor
						for (size_type ll = cst.lb(w); ll <= cst.rb(w); ll++) {
							// get row
							size_type ii = leavesmap[ll];
							assert(leavesmap.count(ll) > 0);
							size_type g = cst.depth(cst.parent(w)) + 1;
							size_type gg = rs_nongaps[ii].rank(x) + g;
							size_type fi;
							if (gg > rs_nongaps[ii].rank(n)) {
								if (!disable_efg_tricks)
									fi = ss_nongaps[ii].select(rs_nongaps[ii].rank(n));
								// fi can be less than x here but it's still correct
								else
									fi = n;
							} else {
								fi = ss_nongaps[ii].select(gg);
							}
							// filter for first occurrence of ignore char
							if (ignorechars and rs_ignores[ii].rank(x) != rs_ignores[ii].rank(n))
								fi = min(ss_ignores[ii].select(rs_ignores[ii].rank(x) + 1), fi);
							if (fi > fimax)
								fimax = fi;
						}
						if (cst.rb(w) == cst.size(cst.root()) - 1)
							break;
						w = cst.select_leaf(cst.rb(w) + 2);
					}
				}
			}
		}
		update_max(f[x], fimax);

		for (size_type i = 0; i < m; i++) {
			for (size_type ll = cst.lb(leaves[i]); ll <= cst.rb(leaves[i]); ll++) { // cst is not a generalized suffix tree
				color[ll] = false;
			}
			if (nongaps[i][x]) {
				leavesmap.erase(cst.lb(leaves[i]));
				leaves[i] = cst.sl(leaves[i]);
				leavesmap[cst.lb(leaves[i])] = i;
				fullrow[i] = false;
			}
		}
	}
}

/* version of compute_f_range that computes range [startx..endx] of f
 * notes: we avoid consuming O(n) bits of space and use O(m) space instead, increasing time? */
template<typename T = size_type>
void compute_f_range(
		size_type const m,
		size_type const n,
		const msa_index &index,
		size_type const startx,
		size_type const endx,
		vector<T> &f, // store result here
		const bool ignorechars,
		const bool disable_efg_tricks
) {
	// unpacking index
	const auto &cst = index.cst;
	const auto &rs_concat_separators = index.rs_concat_separators;
	const auto &nongaps = index.nongaps;
	const auto &rs_nongaps = index.rs_nongaps;
	const auto &ss_nongaps = index.ss_nongaps;
	const auto &ignores = index.ignores;
	const auto &rs_ignores = index.rs_ignores;
	const auto &ss_ignores = index.ss_ignores;

	// find leaves corresponding to suffixes starting at col x+1
	vector<node_t> leaves(m, cst.root());
	unordered_map<size_type, size_type>  leavesmap; // leaf index -> MSA row
	for (size_type indexpos = 0, i = 0; i < m; i++) {
		indexpos += rs_nongaps[i].rank(startx);
		leaves[i] = cst.select_leaf(cst.csa.isa[indexpos] + 1);
		if (rs_nongaps[i].rank(startx) != 0) {
			leavesmap[cst.lb(leaves[i])] = i;
		}
		indexpos += rs_nongaps[i].rank(n) - rs_nongaps[i].rank(startx) + 1;
	}

	for (size_type x = startx; x <= endx; x++) {
#ifdef ALGO_HPP_DEBUG
		cerr << "DEBUG: computing f[x], with x = " << x << endl;
		cerr << "DEBUG: leaves for f[" << x << "] are: ";
		for (auto l : leaves) cerr << l << " ";
		cerr << endl;
#endif
		size_type fimax = x;

		// Process each set of contiguous leaves
		for (size_type i = 0; i < m; i++) {
			node_t const l = leaves[i];
			if (!disable_efg_tricks and rs_nongaps[i].rank(x) == 0) // leaves corresponding to full rows
				continue;
			if (ignorechars and ignores[i][x]) // leaves corresponding to ignorechars
				continue;
			if (cst.lb(l) == 0 || leavesmap.find(cst.lb(l) - 1) == leavesmap.end()) {
				// if leftmost leaf does not correspond to row i, skip
				if (rs_concat_separators.rank(cst.sn(cst.select_leaf(cst.lb(l) + 1))) != i)
					break;
				size_type lb = cst.lb(l);
				size_type rb = cst.rb(l);
				while (rb < cst.size(cst.root()) - 1 && leavesmap.find(cst.rb(cst.select_leaf(rb)) + 2) != leavesmap.end()) {
					rb = rb + 1;
				}

				// Find the exclusive ancestors
				node_t w = l;
				while (cst.rb(w) <= rb) {
					node_t parent = cst.parent(w);
					if (lb <= cst.lb(parent) && cst.rb(parent) <= rb) {
						// parent is a correct replacement
						w = parent;
					} else {
						// parent fails so w is an exclusive ancestor
						for (size_type ll = cst.lb(w); ll <= cst.rb(w); ll++) {
							// get row
							size_type ii = leavesmap[ll];
							assert(leavesmap.count(ll) > 0);
							size_type g = cst.depth(cst.parent(w)) + 1;
							size_type gg = rs_nongaps[ii].rank(x) + g;
							size_type fi;
							if (gg > rs_nongaps[ii].rank(n)) {
								if (!disable_efg_tricks)
									fi = ss_nongaps[ii].select(rs_nongaps[ii].rank(n));
								// fi can be less than x here but it's still correct
								else
									fi = n;
								// fi can be less than x here but it's still correct
							} else {
								fi = ss_nongaps[ii].select(gg);
							}
							// filter for first occurrence of ignore char
							if (ignorechars and rs_ignores[ii].rank(x) != rs_ignores[ii].rank(n))
								fi = min(ss_ignores[ii].select(rs_ignores[ii].rank(x) + 1), fi);
							if (fi > fimax)
								fimax = fi;
						}
						if (cst.rb(w) == cst.size(cst.root()) - 1)
							break;
						w = cst.select_leaf(cst.rb(w) + 2);
					}
				}
			}
		}
		update_max(f[x], fimax);

		for (size_type i = 0; i < m; i++) {
			if (nongaps[i][x]) {
				leavesmap.erase(cst.lb(leaves[i]));
				leaves[i] = cst.sl(leaves[i]);
				leavesmap[cst.lb(leaves[i])] = i;
			}
		}
	}
}

/* version of compute_f that parallelizes the computation of f */
void compute_f_multithread(
		const unsigned long long m,
		const unsigned long long n,
		const msa_index &index,
		const int threads,
		vector<size_type> &f, // store result here
		const bool ignorechars = false,
		const bool disable_efg_tricks = false
) {
	assert(threads > 0);
	const size_type step = n / threads + ((n % threads) > 0);
	unsigned long long consumed = 0;

	vector<thread> t;
	for (int k = 0; k < threads and consumed < n; k++) {
		t.push_back(thread(compute_f_range<>,
				m,
				n,
				ref(index),
				consumed,
				min(n, consumed + step) - 1,
				ref(f),
				ignorechars,
				disable_efg_tricks));

		consumed = min(n, consumed + step);
	}

	for (auto &tt : t) tt.join();
}

/* modifies: destroys f after use
 * notes: Algorithm 1 from 10.1016/j.tcs.2023.114269 */
template<typename T = size_type>
segmentation minmaxlength(
	const unsigned long long n,
	vector<T> &f
) {
#ifdef ALGO_HPP_DEBUG
	cerr << "\nDEBUG: f is = ";
	for (const auto &v : f) cerr << " " << v;
	cerr << endl;
#endif
	// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
	vector<pair<size_type,size_type>> minimal_right_extensions;
	minimal_right_extensions.resize(n);
	for (size_type x = 0; x < n; x++) {
		pair<size_type,size_type> p(x, f[x]+1);
		minimal_right_extensions[x] = p;
	}
	std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), [](std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) { return (std::get<1>(a) < std::get<1>(b)); });
	f.clear();

	// Compute optimal segmentation
	// TODO: swap size_type with int32 or optimal multiple of 2
	vector<size_type> count_solutions(n, 0);
	vector<size_type> backtrack_count(n, 0);
	vector<list<pair<size_type,size_type>>> transition_list(n + 2);
	vector<size_type> minmaxlength(n + 1, 0);
	vector<size_type> backtrack(n + 1, 0);
	size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
	for (size_type j = 1; j <= n; j++) {
		while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
			size_type xy  = std::get<0>(minimal_right_extensions[y]);
			size_type rec_score = minmaxlength[xy];
			if (rec_score > n) {
				// filter out cases when there is no recursive solution
			} else if (j <= xy + rec_score) {
				count_solutions[rec_score] += 1;
				I = std::min(I, rec_score);
				const size_type current_x = backtrack_count[rec_score];
				if (xy + rec_score > current_x + minmaxlength[current_x]) {
					backtrack_count[rec_score] = xy;
				}
				if (xy + rec_score + 1 <= n) {
					transition_list[xy + rec_score + 1].push_back(minimal_right_extensions[y]);
				}
			} else {
				if (j - xy < S) {
					backtrack_S = xy;
				}
				S = std::min(S, j - xy);
			}
			y += 1;
		}
		for (auto pair : transition_list[j]) {
			const size_type x = std::get<0>(pair);
			count_solutions[minmaxlength[x]] -= 1;
			if (j - x < S) {
				S = j - x;
				backtrack_S = x;
			}
			if (count_solutions[minmaxlength[x]] == 0) {
				backtrack_count[minmaxlength[x]] = 0;
			}
		}
		if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S?
			minmaxlength[j] = I;
			backtrack[j] = backtrack_count[I];
		} else {
			minmaxlength[j] = S;
			backtrack[j] = backtrack_S;
		}
		S += 1;
		if (count_solutions[I] == 0)
			I += 1;
	}
#ifdef ALGO_HPP_DEBUG
	cerr << "DEBUG: optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;
#endif

	segmentation SS;
	SS.push_back(n);
	for (size_type j = n; backtrack[j] != 0; j = backtrack[j])
		SS.push_back(backtrack[j]);
	SS.push_back(0);
	std::reverse(SS.begin(), SS.end());

	return SS;
}

/* requires: index to have been generated with function index_external_memory */
segmentation optimal_segmentation_minmaxlength(
		const unsigned long long m,
		const unsigned long long n,
		const msa_index &index,
		const int threads = -1,
		const bool ignorechars = false,
		const bool disable_efg_tricks = false
) {
	//f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is semi-repeat-free
	vector<size_type> f(n, 0);
	if (threads == -1) {
		// single-threaded
		compute_f(m, n, index, f, ignorechars, disable_efg_tricks);
	} else {
		compute_f_multithread(m, n, index, threads, f, ignorechars, disable_efg_tricks);
	}

	if (disable_efg_tricks and f[0] == n) {
		std::cerr << "\n" << "ERROR: no valid segmentation found!" << endl;
		exit(1);
	} // else there is always a valid segmentation

	return minmaxlength(n, f);
}

bool is_valid(
		const efg &g,
		const efg_index &gi,
		const unsigned long node,
		const string &ignorechars
) {
	// unpacking index
	const auto &index = gi.concat;
	const auto &node_labels = gi.labels;
	const auto &node_to_block = g.node_to_block;
	const auto &edge_targets = gi.edge_targets;
	const auto &is_source = gi.is_source;
	const auto &is_sink = gi.is_sink;
	const auto &rs_concat_leaders = gi.rs_concat_leaders;
	const auto &rs_concat_separators = gi.rs_concat_separators;
	const auto &ss_concat_separators = gi.ss_concat_separators;
#ifdef ALGO_HPP_DEBUG
	cerr << "DEBUG: processing node " << node << " : " << node_labels[node] << endl;
#endif

	if (is_source[node] or is_sink[node])
		return true;

	if (ignorechars.length() > 0 and contains_chars(node_labels[node], ignorechars))
		return true;

	sa_type::size_type l, r;
	sdsl::backward_search(index, 0, index.size()-1, node_labels[node].begin(), node_labels[node].end(), l, r);
	assert(l <= r);

	for (sa_type::size_type i = l; i <= r; i++) {
		// locate edge
		const unsigned long occ = index[i];
		const unsigned long occedge = rs_concat_separators(occ);
		const unsigned long occedgeindex = occ - ((occedge == 0) ? 0 : ss_concat_separators(occedge) + 1);
		const unsigned long snode = rs_concat_leaders(occ+1) - 1;
		const unsigned long slength = node_labels[snode].size();

		// locate specific node in the edge
		unsigned long occnode, occnodeindex;
		if (occedgeindex < slength) { // source node
			occnode = snode;
			occnodeindex = occedgeindex;
		} else { // target node
			occnode = edge_targets[occedge];
			occnodeindex = occedgeindex - slength;
		}

#ifdef ALGO_HPP_DEBUG
		cerr << "DEBUG: occurrence is inside edge " << occedge;
		cerr << ", edgeindex " << occedgeindex;
		cerr << ", snode " << snode;
		cerr << ", slength " << slength;
		cerr << ", node " << occnode;
		cerr << ", nodeindex " << occnodeindex;
		cerr << ", block " << node_to_block.at(occnode);
		cerr << endl;
#endif

		// semi-repeat-free property (sources and sinks are special)
		if (occnodeindex != 0) {
			return false;
		}
		if (node_to_block.at(node) != node_to_block.at(occnode)) {
			return false;
		}
		/*if (node != occnode and node_labels[node] == node_labels[occnode]) {
			return false;
		}*/
	}
	return true;
}

void is_block_valid_worker(
		const efg &g,
		const efg_index &gi,
		const string &ignorechars,
		const size_type start_block,
		const int step_block,
		vector<bool> &to_remove // fill in this
) {
	for (size_type b = start_block; b < g.blocks.size(); b += step_block) {
		for (auto &[_, node] : g.blocks[b]) {
			if (!is_valid(g, gi, node, ignorechars)) {
				to_remove[b] = true;
				continue;
			}
		}
	}
}

/* requires: graph g to be built from segmentation S with function segment_msa
 * modifies: g is transformed into graph/index pair (g, gi), see index_efg_external_memory
 *           to_remove is replaced and marks the invalid segments
 * returns: true if all segments are semi-repeat-free
 * note: sources and sinks are implicitly unique (see tricks from 10.1093/bioinformatics/btaf225) */
bool is_indexable(
		const segmentation &S,
		efg &g,
		const string &ignorechars,
		const int threads,
		const path &tmpdir,
		vector<bool> &to_remove
) {
	bool res = true;
	efg_index gi = index_efg_external_memory(g, tmpdir);
	to_remove = vector<bool>(g.blocks.size(), false);
	if (threads == -1) {
		for (size_type b = 0; b < g.blocks.size(); b++) {
			for (auto &[_, node] : g.blocks[b]) {
				if (!is_valid(g, gi, node, ignorechars)) {
					to_remove[b] = true;
					res = false;
					continue;
				}
			}
		}
	} else {
		assert(threads >= 1);
		vector<thread> t;
		for (int k = 0; k < threads; k++) {
			t.push_back(thread(is_block_valid_worker,
						ref(g),
						ref(gi),
						ref(ignorechars),
						k,
						threads,
						ref(to_remove)));
		}
		for (auto &tt : t) tt.join();
		for (auto r : to_remove) {
			if (r) res = false;
		}
	}
#ifdef ALGO_HPP_DEBUG
	cerr << "DEBUG: to_remove is";
	for (const auto bb : to_remove) cerr << " " << ((bb) ? 1 : 0);
	cerr << endl;
#endif
	return res;
}

void prune(segmentation &S, vector<bool> &to_remove)
{
	// S is 1 unit longer than to_remove
	assert(S.size() == to_remove.size() + 1);

	segmentation S_new;
	S_new.push_back(S.front());
	for (size_type i = 1; i < S.size() - 1; i++) {
		if (!to_remove[i])
			S_new.push_back(S[i]);
	}
	S_new.push_back(S.back());
	std::swap(S_new, S);
}

const string _ignorechars = "";
void heuristic_index_compute_f_worker(
		const path &tmpdir,
		const int thread_id,
		const int heuristic_subset,
		atomic<unsigned long long> &m, // update this to total rows
		unsigned long long &n,
		vector<atomic<size_type>> &f,
		const string &ignorechars = _ignorechars,
		const bool disable_efg_tricks = false
) {
	const path &thread_tmpdir = tmpdir / ("thread_" + std::to_string(thread_id+1));
	if (exists(thread_tmpdir)) {
		cerr << "ERROR (thread " << thread_id+1 << "): directory " << thread_tmpdir << " already exists! Please clean up the temporary directory." << endl;
		exit(1);
	}
	create_directory(thread_tmpdir);
	unsigned long long mm = 0;
	unsigned long long startrow = 0;
	do {
		msa_index index = index_external_memory(path(), thread_tmpdir, ignorechars, mm, n, startrow, heuristic_subset);
		if (mm == 0) break;
		compute_f_range<atomic<size_type>>(mm, n, index, 0, n - 1, f, (ignorechars != ""), disable_efg_tricks);
		// alternatively, compute_f<atomic<size_type>>(mm, n, index, f, (ignorechars != ""), disable_efg_tricks);
#ifdef ALGO_HPP_DEBUG
		double fmean = 0;
		for (const auto &ff : f) fmean += ff;
		fmean /= f.size();
		cerr << "\nDEBUG: fmean is " << fmean << endl;
#endif
		m += mm;
	} while (true);
	remove(thread_tmpdir);
}

/* notes: streams from disk heuristic_subset (times the number of threads) rows at a time */
segmentation heuristic_subset_segmentation_minmaxlength(
		const path &msa_path,
		const path &tmpdir,
		const int heuristic_subset,
		unsigned long long &m, // update this
		unsigned long long &n, // update this
		const int threads = -1,
		const string &ignorechars = _ignorechars,
		const bool disable_efg_tricks = false
) {
	assert(heuristic_subset > 0);

	string line, _;
	{
		open_msa_file(msa_path);
		if (!get_msa_line(line, _)) { cerr << "ERROR: cannot read MSA line!" << endl; exit(1); };
		n = line.size();
		close_msa_file();
	}

	if (threads == -1) {
		open_msa_file(msa_path);
		// single-threaded
		// f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is semi-repeat-free
		vector<size_type> f(n, 0);
		unsigned long long mm = 0;
		unsigned long long startrow = 0;
		do {
			msa_index index = index_external_memory(path(), tmpdir, ignorechars, mm, n, startrow, heuristic_subset);
			if (mm == 0) break;
			compute_f(mm, n, index, f, (ignorechars != ""), disable_efg_tricks);
#ifdef ALGO_HPP_DEBUG
			double fmean = 0;
			for (const auto &ff : f) fmean += ff;
			fmean /= f.size();
			cerr << "\nDEBUG: fmean is " << fmean << endl;
#endif
			startrow += mm;
		} while (true);
		close_msa_file();
		return minmaxlength(n, f);
	} else {
		open_msa_file(msa_path);
		// multi-threaded
		// f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is semi-repeat-free
		vector<atomic<size_type>> f(n);
		for (size_type i = 0; i < f.size(); i++) f[i] = 0;
		atomic<unsigned long long> mm = 0;
		vector<thread> t;
		for (int k = 0; k < threads; k++) {
			t.push_back(thread(heuristic_index_compute_f_worker,
						ref(tmpdir),
						k,
						heuristic_subset,
						ref(mm),
						ref(n),
						ref(f),
						ref(ignorechars),
						disable_efg_tricks));
		}
		for (auto &tt : t) tt.join();
		close_msa_file();
		m = mm;
		return minmaxlength(n, f);
	}
}

} // namespace fbg::algo
#endif // ifndef ALGO_HPP
