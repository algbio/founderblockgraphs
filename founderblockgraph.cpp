/**
 * founderblockgraph
 * Copyright (C) 2020-2021 Veli Mäkinen, Tuukka Norri
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <chrono>
#include <iostream>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/util.hpp>
#include <stdexcept>
#include <unistd.h> // mkstemp
#include <unordered_set>
#include "founderblockgraph_cmdline.h"
#include "founder_block_index.hpp"


namespace {
	
	namespace chrono = std::chrono;
	namespace fbg = founder_block_graph;

	typedef sdsl::cst_sct3 <>				cst_type;
	typedef cst_type::size_type				size_type;

	typedef std::unordered_set <size_type>	edge_set;
	typedef std::vector <edge_set>			adjacency_list;
	
	
	inline constexpr size_type const INVALID_SIZE{std::numeric_limits <size_type>::max()};
	
	
	class temporary_file
	{
	protected:
		std::string	m_name;
		int			m_fd{-1};
	
	public:
		temporary_file(std::string const &dst_dir):
			m_name(dst_dir + "/temp.XXXXXX")
		{
		}
	
		~temporary_file()
		{
			if (-1 != m_fd)
				close(m_fd);
		}
		
		std::string const &name() const { return m_name; }
		
		void open()
		{
			m_fd = mkstemp(m_name.data()); // Since C++11, data() returns a null-terminated buffer.
			if (-1 == m_fd)
			{
				std::string reason("Unable to create temporary file: ");
				reason += strerror(errno);
				throw std::runtime_error(reason);
			}
		}
	};
	
	
	template <typename t_string>
	std::string remove_gaps(t_string const &S) {
	    std::string result = "";
	    for (size_type i=0; i<S.size(); i++) 
	        if (S[i]!='-')
	            result += S[i];
	    return result;
	}
	
	
	bool check_gaps(std::string const &identifier, std::string const &sequence, std::size_t gap_limit)
	{
		if (0 == gap_limit)
			return true;
		
		std::size_t gaprun(0);
		std::size_t maxgaprun(0);
		for (auto const c : sequence)
		{
			if ('-' == c)
				++gaprun;
			else
			{
				maxgaprun = std::max(gaprun, maxgaprun);
				gaprun = 0;
			}
		}
		
		maxgaprun = std::max(gaprun, maxgaprun);
		if (maxgaprun < gap_limit)
			return true;
		else
		{
			std::cerr << "NOTICE: Sequence “"
				<< identifier.substr(1)
				<< "” contained a gap run with "
				<< maxgaprun
				<< " characters.\n";
			return false;
		}
	}
	
	
	bool check_sequence_length(std::string const &identifier, std::string const &seq, std::size_t expected_length)
	{
		if (seq.size() == expected_length)
			return true;
		
		std::cerr << "WARNING: length of the sequence “"
			<< identifier.substr(1)
			<< "” does not match that of the first sequence; skipping. ("
			<< expected_length
			<< " vs. "
			<< seq.size()
			<< ")\n";
		return false;
	}
	
	
	void read_input(char const *input_path, std::size_t gap_limit, std::vector <std::string> &msa)
	{
		std::string line, identifier, entry;

		// Reading input fasta
		std::fstream fs;
		fs.open(input_path, std::fstream::in);
		
		// Assume that the first line contains a header.
		std::getline(fs, identifier);
		
		std::size_t expected_length(0);
		bool is_first(true);
		while (std::getline(fs, line))
		{
			if (line[0] == '>') // header
			{
				if (is_first)
				{
					expected_length = entry.size();
					is_first = false;
				}
				
				if (check_sequence_length(identifier, entry, expected_length) &&
					check_gaps(identifier, entry, gap_limit))
				{
					msa.push_back(entry);
				}
				entry.clear();
				identifier = line;
			}
			else
			{
				entry += line;
			}
		}
		fs.close();
		if (check_sequence_length(identifier, entry, expected_length) &&
			check_gaps(identifier, entry, gap_limit))
		{
			msa.push_back(entry);
		}
	}
	
	
	bool load_cst(
		char const *input_path_c,
		std::vector <std::string> const &msa,
		cst_type &reverse_cst,
		std::size_t gap_limit
	)
	{
		std::string const input_path(input_path_c);
		std::string const index_suffix(".cst");
		std::string const plain_suffix(".plain");
		std::string const reverse_suffix(".reverse");
		auto const gap_limit_s(std::to_string(gap_limit));
		
		auto const reverse_index_file(input_path + reverse_suffix + plain_suffix + gap_limit_s + index_suffix);
		
		// Construct compressed suffix trees for C.
		if (!sdsl::load_from_file(reverse_cst, reverse_index_file))
		{
			std::cout << "No index " << reverse_index_file << " located. Building index now.\n";
			auto const concatenated_inputs_path(input_path + reverse_suffix + plain_suffix);
			
			// Output the exact reverse of the above to disk.
			{
				std::fstream fs;
				fs.open(concatenated_inputs_path, std::fstream::out);
				auto seq_it(msa.rbegin());
				auto const seq_end(msa.rend());
				while (seq_it != seq_end)
				{
					fs << '#';
					
					auto it(seq_it->rbegin());
					auto const end(seq_it->rend());
					while (it != end)
					{
						auto const c(*it);
						if ('-' != c)
							fs << c;
						++it;
					}
					++seq_it;
				}
				
				fs.close();
			}
			
			std::ifstream in(concatenated_inputs_path);
			if (!in)
			{
				std::cout << "ERROR: File " << concatenated_inputs_path << " does not exist.\n";
				return false;
			}
			sdsl::construct(reverse_cst, concatenated_inputs_path, 1); // generate index
			sdsl::store_to_file(reverse_cst, reverse_index_file); // save it
		}

		return true;
	}
	
	
	struct interval
	{
		size_type sp{};
		size_type ep{};

		interval() = default;
		
		interval(size_type l, size_type r) : sp(l), ep(r) {}
		
		cst_type::node_type find_parent(cst_type const &cst);
		size_type length() const { return ep - sp + 1; }
		inline bool operator<(interval const &other) const;
	};
	
	
	// Support for contract_right() by one character.
	class interval_pair
	{
	public:
		typedef cst_type::node_type node_type;
		
	protected:
		interval	m_parent;
		interval	m_current;
		size_type	m_edge_characters_remaining{};
		size_type	m_suffix_length{};
		
	public:
		interval_pair() = default;
		
		interval_pair(cst_type const &rev_cst):
			m_parent(INVALID_SIZE, 0),
			m_current(0, rev_cst.csa.size() - 1)
		{
		}
		
		size_type sp() const { return m_current.sp; }
		size_type ep() const { return m_current.ep; }
		size_type length() const { return m_current.length(); }
		interval get_interval() const { return m_current; }
		size_type get_suffix_length() const { return m_suffix_length; }
		void extend_left(cst_type const &cst, char const cc);
		void contract_right(cst_type const &cst);
	};
	
	
	bool interval::operator<(interval const &other) const
	{
		if (sp < other.sp) return 1;
		else if (sp == other.sp && ep > other.ep) return 1; // *this contains other but is not equivalent with it.
		else return 0;
	}
	
	
	auto interval::find_parent(cst_type const &cst) -> cst_type::node_type
	{
		auto const ll(cst.select_leaf(sp + 1));
		auto const rl(cst.select_leaf(ep + 1));
		auto const w(cst.lca(ll, rl));
		auto const u(cst.parent(w));
		return u;
	}


	void extract(char const *msg, cst_type const &cst, cst_type::node_type const &node)
	{
		std::cerr << msg;
		auto const str(sdsl::extract(cst, node));
		for (auto const c : str)
			std::cerr << char(c);
		std::cerr << '\n';
	}


	void extract(char const *msg, cst_type const &cst, interval const &interval, size_type const length)
	{
		auto const ll(cst.select_leaf(interval.sp + 1));
		auto const rl(cst.select_leaf(interval.ep + 1));
		auto const node(cst.lca(ll, rl));
		
		std::cerr << msg;
		if (length)
		{
			auto const str(sdsl::extract(cst, node));
			for (size_type i(0); i < length; ++i)
				std::cerr << char(str[i]);
		}
		std::cerr << '\n';
	}


	void interval_pair::extend_left(cst_type const &cst, char const cc)
	{
		//std::cerr << "------\n";

		// Mark the parent invalid.
		m_parent.sp = INVALID_SIZE;
		
		// Update the lexicographic range.
		auto &sp(m_current.sp);
		auto &ep(m_current.ep);
		sdsl::backward_search(cst.csa, sp, ep, cc, sp, ep);

		++m_suffix_length;
		
		//extract("node after extend_left:       ", cst, m_current, m_suffix_length);
	}
	
	
	void interval_pair::contract_right(cst_type const &cst)
	{
		//std::cerr << "------\n";
		
		// While the represented suffix corresponds to a position
		// in the middle of a CST edge, the lexicographic range remains
		// the same.
		
		assert(m_suffix_length); // Can only right-contract non-root.
		
		if (INVALID_SIZE == m_parent.sp)
		{
			// Parent is currently not valid; find it.
			auto parent_node(m_current.find_parent(cst));
			size_type parent_depth(0);
			while (true)
			{
				parent_depth = cst.depth(parent_node);
				if (parent_depth < m_suffix_length)
					break;
				
				// Make sure the loop is finite.
				if (parent_node == cst.root())
					throw std::logic_error("Reached the root node but still trying to find a suitable parent");
				
				assert(parent_node != cst.root());
				parent_node = cst.parent(parent_node);
			}
			m_parent.sp = cst.lb(parent_node);
			m_parent.ep = cst.rb(parent_node);
			m_edge_characters_remaining = m_suffix_length - parent_depth;
		}
		
		// Check whether the parent was reached.
		assert(m_edge_characters_remaining);
		--m_edge_characters_remaining;
		--m_suffix_length;
		if (0 == m_edge_characters_remaining)
		{
			m_current = m_parent;
			m_parent.sp = INVALID_SIZE;
		}
		//extract("node after contract_right:    ", cst, m_current, m_suffix_length);
	}
	
	typedef std::vector <interval> interval_vector;
	typedef std::vector <interval_pair> interval_pair_vector;
	
	
	struct column_range
	{
		size_type lb{};
		size_type rb{};
		
		column_range() = default;
		
		column_range(size_type lb_, size_type rb_):
			lb(lb_),
			rb(rb_)
		{
		}
		
		size_type length() const { return rb - lb + 1; }
		size_type extended_length(size_type lb_) const { return rb - lb_ + 1; }
	};
	
	typedef std::vector <column_range> column_range_vector;
	
	
	size_type sum_interval_lengths(interval_vector const &sorted_bwt_intervals)
	{
		size_type retval(0);
		auto it(sorted_bwt_intervals.begin());
		auto const end(sorted_bwt_intervals.end());
		while (it != end)
		{
			retval += it->length();
			auto const ep(it->ep);
			
			// Find the first interval not enclosed by *it.
			it = std::partition_point(it, end, [ep](interval const &inv) {
				return inv.ep <= ep;
			});
		}
		
		return retval;
	}
	
	
	void extend_ranges_left(
		std::vector <std::string> const &msa,
		cst_type const &cst,
		size_type const input_pos,
		interval_pair_vector &bwt_intervals
	)
	{
		size_type const m(msa.size());
		
		// Find the next character on each input row.
		for (size_type i(0); i < m; ++i)
		{
			auto const cc(msa[i][input_pos]);
			if ('-' != cc)
				bwt_intervals[i].extend_left(cst, cc);
		}
	}


	void contract_ranges_right(
		std::vector <std::string> const &msa,
		cst_type const &cst,
		size_type const input_pos,
		interval_pair_vector &bwt_intervals
	)
	{
		size_type const m(msa.size());

		// For each interval, check whether it can be contracted.
		for (size_type i(0); i < m; ++i)
		{
			auto const cc(msa[i][input_pos]);
			if ('-' == cc)
				continue;
			
			// The current character was not a gap; contract.
			bwt_intervals[i].contract_right(cst);
		}
	}
	
	
	bool check_block(
		size_type const m,
		interval_pair_vector const &bwt_intervals,
		interval_vector &sorted_bwt_intervals
	)
	{
		// Copy the intervals for sorting.
		// (We could maintain just a permutation but this is easier.)
		sorted_bwt_intervals.clear();
		std::transform(bwt_intervals.begin(), bwt_intervals.end(), std::back_inserter(sorted_bwt_intervals),
					   [](interval_pair const &inv){
			return inv.get_interval();
		});
		
		std::sort(sorted_bwt_intervals.begin(), sorted_bwt_intervals.end());
		auto const length_sum(sum_interval_lengths(sorted_bwt_intervals));
		
		assert(m <= length_sum);
		return (m == length_sum);
	}


	void find_valid_blocks(
		std::vector <std::string> const &msa,
		cst_type const &rev_cst,
		column_range_vector &column_ranges
	)
	{
		size_type const m(msa.size());
		assert(0 < m);
		size_type const n(msa[0].size());
		
		// Initially all of the intervals correspond to all suffixes of the input.
		interval_pair_vector bwt_intervals(m, interval_pair(rev_cst));
		interval_pair_vector prev_bwt_intervals;
		interval_vector sorted_bwt_intervals;
		
		size_type column_range_lb(0);
		size_type column_range_rb(0); // Closed interval.
		bool found_valid_range(false);
		
		// Process the input from right to left.
		// We extend the window right and contract it from the left, since
		// the index is for the reverse of the text (while the names of the
		// corresponding functions refer to the original text).
		// First try to find a suitable initial range.
		while (column_range_rb < n)
		{
			extend_ranges_left(msa, rev_cst, column_range_rb, bwt_intervals);
			if (check_block(m, bwt_intervals, sorted_bwt_intervals))
			{
				column_ranges.emplace_back(column_range_lb, column_range_rb);
				found_valid_range = true;
				break;
			}
			
			++column_range_rb;
		}
		
		if (!found_valid_range)
			return;

		// Continue by right-contracting the initial column range to find the shortest
		// range that is segment repeat-free. Then left-extend by one and repeat.
		while (true)
		{
			while (true)
			{
				// Try to left-contract the BWT intervals.
				prev_bwt_intervals = bwt_intervals;
				assert(column_range_lb <= column_range_rb);
				contract_ranges_right(msa, rev_cst, column_range_lb, bwt_intervals);
				
				if (check_block(m, bwt_intervals, sorted_bwt_intervals))
				{
					// The new interval was valid; store it.
					// The previous value may be overwritten since the new range is shorter.
					++column_range_lb;
					auto &last(column_ranges.back());
					if (last.rb == column_range_rb)
						last.lb = column_range_lb;
					else
						column_ranges.emplace_back(column_range_lb, column_range_rb);
				}
				else
				{
					// Restore the previous state.
					bwt_intervals = prev_bwt_intervals;
					break;
				}
			}
			
			++column_range_rb;
			if (n <= column_range_rb)
				break;
			
			// Right-extend.
			extend_ranges_left(msa, rev_cst, column_range_rb, bwt_intervals);
		}
	}
	
	
	void segment(
		std::vector <std::string> const &msa,
		cst_type const &rev_cst,
		std::vector <std::string> &out_labels,
		adjacency_list &out_edges
	)
	{
		assert(!msa.empty());
		size_type const n(msa[0].size());
		size_type const m(msa.size());
		
		// Repeat-free segments.
		std::vector <column_range> column_ranges;
		find_valid_blocks(msa, rev_cst, column_ranges);
		
		auto const &rightmost_range(column_ranges.back());
		if (rightmost_range.rb != n - 1)
		{
			std::cerr << "Only the trivial segmentation is valid.\n";
			return;
		}

		// Verify the monotonicity of the column ranges.
		assert(std::is_sorted(column_ranges.begin(), column_ranges.end(), [](column_range const &lhs, column_range const &rhs){
			return lhs.lb <= rhs.lb && lhs.rb < rhs.rb;
		}));
		
		auto const range_count(column_ranges.size());
		std::vector <size_type> s(n, INVALID_SIZE);
		std::vector <size_type> next(range_count, INVALID_SIZE); // Next non-overlapping block indices
		std::vector <size_type> trace(range_count, INVALID_SIZE); // Previous block indices in the optimal path
				
		// Initial scores.
		for (auto const &range : column_ranges)
		{
			assert(range.rb < s.size());
			s[range.rb] = range.extended_length(0);
		}
		
		// Next ranges that do not overlap with the current one.
		{
			size_type i(0);
			size_type j(0);
			
			while (i < range_count)
			{
				assert(i < column_ranges.size());
				auto const &lhs(column_ranges[i]);
				while (j < range_count)
				{
					assert(j < column_ranges.size());
					auto const &rhs(column_ranges[j]);
					if (lhs.rb < rhs.lb)
					{
						assert(i < next.size());
						next[i] = j;
						break;
					}
					++j;
				}
				++i;
			}
		}
		
		// Determine the score limit.
		size_type score_limit(0);
		{
			auto const &leftmost_range(column_ranges.front());
			size_type max_score(leftmost_range.extended_length(0));
			
			for (size_type i(0); i < range_count; ++i)
			{
				assert(i < column_ranges.size());
				assert(i < next.size());
				auto const &lhs_range(column_ranges[i]);
				max_score = std::max(max_score, lhs_range.length());
				
				auto const next_idx(next[i]);
				if (INVALID_SIZE == next_idx)
					continue;
				
				assert(next_idx < column_ranges.size());
				auto const &rhs_range(column_ranges[next_idx]);
				max_score = std::max(max_score, rhs_range.extended_length(lhs_range.rb + 1)); // Allow extensions.
			}
			score_limit = max_score + 1;
		}
		
		// Main algorithm.
		{
			size_type i(0);
			auto const range_count_1(range_count - 1);
			while (i < range_count_1)
			{
				assert(i < column_ranges.size());
				auto const &lhs(column_ranges[i]);
				
				// Compare to the ranges to the right starting from next[i].
				assert(i < next.size());
				for (auto j(next[i]); j < range_count; ++j)
				{
					assert(j < column_ranges.size());
					auto const &rhs(column_ranges[j]);
					assert(lhs.rb < rhs.lb);
					auto const rhs_score(rhs.extended_length(lhs.rb + 1));
					
					// Stop if the score limit has been reached.
					if (score_limit <= rhs_score)
						break;
					
					assert(lhs.rb < s.size());
					auto const score(std::max(s[lhs.rb], rhs_score));
					if (score < s[rhs.rb])
					{
						assert(rhs.rb < s.size());
						assert(j < trace.size());
						s[rhs.rb] = score;
						trace[j] = i;
					}
				}
				
				++i;
			}
		}
		
		// Extend the column ranges on the found path.
		std::vector <column_range> chosen_ranges;
		{
			auto i(range_count - 1);
			while (true)
			{
				assert(i < column_ranges.size());
				assert(i < trace.size());
				auto const prev_idx(trace[i]);
				auto &range(column_ranges[i]);
				if (INVALID_SIZE == prev_idx)
				{
					range.lb = 0;
					chosen_ranges.push_back(range);
					break;
				}
				
				assert(prev_idx < column_ranges.size());
				auto const &prev_range(column_ranges[prev_idx]);
				range.lb = prev_range.rb + 1;
				chosen_ranges.push_back(range);
				i = prev_idx;
			}
			
			std::reverse(chosen_ranges.begin(), chosen_ranges.end());
		}
		
		// Create labels and edges for each column range.
		{
			std::set <std::string> lhs_labels;
			std::set <std::string> rhs_labels;
			
			// Create string_views from the input.
			std::vector <std::string_view> msa_views(m);
			for (size_type i(0); i < m; ++i)
				msa_views[i] = msa[i];
			
			// Create labels for the first range.
			auto const &first_range(chosen_ranges.front());
			for (size_type i(0); i < m; ++i)
			{
				auto const sv(msa_views[i].substr(first_range.lb, first_range.length()));
				lhs_labels.emplace(remove_gaps(sv));
			}
			
			// Copy to the output list.
			std::copy(lhs_labels.begin(), lhs_labels.end(), std::back_inserter(out_labels));
			
			// Create labels and edges for the remaining ranges.
			size_type lhs_idx_start(0);
			size_type rhs_idx_start(lhs_labels.size());
			auto const range_count(chosen_ranges.size());
			for (size_type j(1); j < range_count; ++j)
			{
				// Create the labels.
				assert(j < chosen_ranges.size());
				auto const &rhs_range(chosen_ranges[j]);
				for (size_type i(0); i < m; ++i)
				{
					assert(i < msa_views.size());
					auto const sv(msa_views[i].substr(rhs_range.lb, rhs_range.length()));
					rhs_labels.emplace(remove_gaps(sv));
				}
				
				auto const lc(lhs_labels.size());
				auto const rc(rhs_labels.size());
				
				// Copy to the output list.
				assert(out_labels.size() == lhs_idx_start + lc);
				assert(out_labels.size() == rhs_idx_start);
				std::copy(rhs_labels.begin(), rhs_labels.end(), std::back_inserter(out_labels));
				assert(out_labels.size() == rhs_idx_start + rc);
				
				// Create the edges.
				assert(out_edges.size() == lhs_idx_start);
				for (size_type ll(0); ll < lc; ++ll)
				{
					auto &adj_list(out_edges.emplace_back());
					for (size_type rl(0); rl < rc; ++rl)
						adj_list.insert(rhs_idx_start + rl);
				}
				assert(out_edges.size() == lhs_idx_start + lc);
				
				// Prepare for the next iteration.
				lhs_idx_start = rhs_idx_start;
				rhs_idx_start += rhs_labels.size();
				
				using std::swap;
				swap(lhs_labels, rhs_labels);
				rhs_labels.clear();
			}
			
			// Make sure that there is an adjacency list for every node.
			out_edges.resize(out_labels.size());
		}
	}
	
	
	void make_index(
	    std::vector <std::string> node_labels,
	    adjacency_list edges,
	    char const *dst_path_c,
	    char const *mem_chart_path_c
	) {
	    // Create a temporary file.
	    std::string const dst_path(dst_path_c);
	    std::string const dirname(sdsl::util::dirname(dst_path));
	    temporary_file temp_file(dirname);
	    temp_file.open();
	    std::fstream temp_os;
	    temp_os.exceptions(std::fstream::failbit);
	    temp_os.open(temp_file.name(), std::ios_base::out);
    
	    // Write the index contents to it.
	    // Also keep track of node occurrences.
	    assert(edges.size() == node_labels.size());
	    std::vector <std::size_t> dst_occurrences(node_labels.size(), 0);
	    for (std::size_t i(0), count(node_labels.size()); i < count; ++i)
	    {
	        auto const &src_label(node_labels[i]);
        
	        // Make sure that the destination nodes are sorted.
	        // This is not necessary for the algorithm to work but is useful for debugging purposes.
	        auto const &dst_nodes(edges[i]);
	        std::vector <size_type> sorted_dst_nodes(dst_nodes.begin(), dst_nodes.end());
	        std::sort(sorted_dst_nodes.begin(), sorted_dst_nodes.end());
        
	        for (auto const dst_node : sorted_dst_nodes)
	        {
	            auto const &dst_label(node_labels[dst_node]);
	            temp_os << fbg::g_separator_character << src_label << dst_label;
	            ++dst_occurrences[dst_node];
	        }
	    }
    
	    temp_os << std::flush;
	    temp_os.close();
    
	    // Construct the CSA.
	    fbg::csa_type csa;
	    sdsl::construct(csa, temp_file.name(), 1);
    
	    // Prepare the B and E bit vectors.
	    auto const csa_size(csa.size());
	    sdsl::bit_vector b_positions(csa_size, 0);
	    sdsl::bit_vector e_positions(csa_size, 0);
	    for (std::size_t i(0), count(node_labels.size()); i < count; ++i)
	    {
	        auto const &label(node_labels[i]);
        
	        fbg::csa_type::size_type lhs{}, rhs{};
	        auto const match_count(
	            sdsl::backward_search(
	                csa,
	                0,
	                csa_size - 1,
	                label.begin(),
	                label.end(),
	                lhs,
	                rhs
	            )
	        );

			assert(lhs < b_positions.size());
			assert(rhs < e_positions.size());
        
	        // Check by searching the separator character.
			// FIXME: none of these are useful with the elastic version, hence the use of false in the beginning.
	        if (false && (b_positions[lhs] || e_positions[rhs] || match_count != (edges[i].size() + dst_occurrences[i])))
	        {
				if (b_positions[lhs])
					std::cerr << "b_positions[lhs] was already set.\n";
				if (e_positions[rhs])
					std::cerr << "e_positions[rhs] was already set.\n";
				if (match_count != (edges[i].size() + dst_occurrences[i]))
					std::cerr << "Got unexpected match count.\n";
	            std::cerr << "Match count:         " << match_count << '\n';
				std::cerr << "Node:                " << i << '\n';
				std::cerr << "Label:               “" << label << "”\n";
	            std::cerr << "Actual in-edges:     " << dst_occurrences[i] << '\n';
				std::cerr << "Out-edges:           " << edges[i].size() << "\n";
	            std::cerr << "Lexicographic range: " << lhs << ", " << rhs << '\n';

	            std::vector <char> buffer;
	            std::cerr << "Suffixes:\n";
	            for (size_type i(lhs); i <= rhs; ++i)
	            {
	                auto const text_pos(csa[i]);
	                auto const end_pos(std::min(text_pos + 100, csa.size() - 1));
                
	                buffer.clear();
	                buffer.resize(100, '\0');
	                auto const count(sdsl::extract(csa, text_pos, end_pos, buffer.begin()));
	                std::cerr << "i: " << i << "\ttext: ";
	                for (size_type j(0); j < count; ++j)
	                    std::cerr << buffer[j];
	                std::cerr << '\n';
	            }
	            std::cerr << "Out-edges:\n";
	            for (auto const dst_node : edges[i])
	                std::cerr << dst_node << '\t' << node_labels[dst_node] << '\n';
            
	            {
	                std::cerr << "Preceding characters:\n";
	                auto const sigma(csa.wavelet_tree.sigma);
	                std::vector<fbg::csa_type::wavelet_tree_type::value_type> cs(sigma, 0);
	                size_type symbol_count(0);
	                std::vector <size_type> r1(sigma, 0), r2(sigma, 0);
	                sdsl::interval_symbols(csa.wavelet_tree, lhs, rhs, symbol_count, cs, r1, r2);
	                for (std::size_t j(0); j < symbol_count; ++j)
	                    std::cerr << "‘" << cs[j] << "’ (" << int(cs[j]) << ")\n";
	            }
            
	            // Output a sample of the node labels, too.
	            // We use binary search here instead of constructing a BWT index of the node labels only.
	            {
	                typedef std::pair<std::size_t, std::string> label_pair;
	                typedef std::vector<label_pair> label_pair_vector;
	                label_pair_vector sorted_node_labels(node_labels.size());
	                for (std::size_t j(0), count(node_labels.size()); j < count; ++j)
	                {
	                    sorted_node_labels[j] = label_pair(j, node_labels[j]);
	                    std::sort(sorted_node_labels.begin(), sorted_node_labels.end(), [](auto const &lhs, auto const &rhs){
	                        return lhs.second < rhs.second;
	                    });
	                }
                
	                // Find the problematic label.
	                auto const range(std::equal_range(sorted_node_labels.begin(), sorted_node_labels.end(), label_pair(0, label), [](auto const &lhs, auto const &rhs){
	                    return lhs.second < rhs.second;
	                }));
	                typedef label_pair_vector::difference_type diff_type;
	                auto it(range.first - std::min(diff_type(5), std::distance(sorted_node_labels.begin(), range.first)));
	                auto const end(range.second + std::min(diff_type(5), std::distance(range.second, sorted_node_labels.end())));
                                 
	                std::cerr << "Sample of node labels:\n";
	                while (it != end)
	                {
	                    std::cerr << it->first << '\t' << it->second << '\n';
	                    ++it;
	                }
	            }
            
	            abort();
	        }
        
	        b_positions[lhs] = 1;
	        e_positions[rhs] = 1;
	    }
    
	    fbg::founder_block_index founder_block_index(
	        std::move(csa),
	        std::move(b_positions),
	        std::move(e_positions)
	    );
    
	    // Write the memory chart if needed.
	    if (mem_chart_path_c)
	    {
	        std::fstream mem_chart_os;
	        mem_chart_os.exceptions(std::fstream::failbit);
	        mem_chart_os.open(mem_chart_path_c, std::ios_base::out);
	        sdsl::write_structure<sdsl::HTML_FORMAT>(founder_block_index, mem_chart_os);
	        mem_chart_os.close();
	    }
    
	    // Write the index.
	    std::fstream index_os;
	    index_os.exceptions(std::fstream::failbit);
	    index_os.open(dst_path, std::ios_base::out);
	    founder_block_index.serialize(index_os);
	    index_os.close();
	}


	void output_graphviz_label(std::ostream &os, std::string const &label) {
	    for (auto const c : label)
	    {
	        if ('"' == c)
	            os << "\\\"";
	        else
	            os << c;
	    }
	}


	void output_graphviz(
	    std::vector <std::string> node_labels,
	    adjacency_list edges,
	    char const *dst_path
	) {
	    std::fstream os;
	    os.exceptions(std::fstream::failbit);
	    os.open(dst_path, std::ios_base::out);
    
	    os << "digraph founder_block_graph {\n";
	    os << "rankdir=\"LR\"\n";
    
	    // Node labels.
	    {
	        std::size_t i(0);
	        for (auto const &label : node_labels)
	        {
	            os << 'n' << i << " [label = \"";
	            output_graphviz_label(os, label);
	            os << "\"];\n";
	            ++i;
	        }
	    }
    
	    // Edges.
	    {
	        auto const count(edges.size());
	        for (std::size_t i(0); i < count; ++i)
	        {
	            auto const &dst_nodes(edges[i]);
	            if (!dst_nodes.empty())
	            {
	                os << 'n' << i << " -> {";
	                bool is_first(true);
	                for (auto const dst_node : dst_nodes)
	                {
	                    if (is_first)
	                        is_first = false;
	                    else
	                        os << " ; ";
	                    os << 'n' << dst_node;
	                }
	                os << "}\n";
	            }
	        }
	    }
    
	    os << "}\n";
	    os.close();
	}
}


int main(int argc, char **argv) { 

	gengetopt_args_info args_info;
	if (0 != cmdline_parser(argc, argv, &args_info))
		return EXIT_FAILURE;
	
	std::ios_base::sync_with_stdio(false);    // Don't use C style IO after calling cmdline_parser.
	
	auto const gap_limit(args_info.gap_limit_arg);
	if (gap_limit < 0)
	{
		std::cerr << "Gap limit needs to be non-negative.\n";
		return EXIT_FAILURE;
	}
	
	auto const start(chrono::high_resolution_clock::now());
	
	std::vector <std::string> msa;
	read_input(args_info.input_arg, gap_limit, msa);
	if (msa.empty())
	{
		std::cerr << "Unable to read sequences from the input\n.";
		return EXIT_FAILURE;
	}
	std::cout << "Input MSA[1.." << msa.size() << ",1.." << msa[0].size() << "]\n";
	
	cst_type reverse_cst;
	if (!load_cst(args_info.input_arg, msa, reverse_cst, gap_limit))
		return EXIT_FAILURE;
	
	std::cout << "Index construction complete, index requires " << sdsl::size_in_mega_bytes(reverse_cst) << " MiB.\n";
	
	std::vector <std::string> node_labels;
	adjacency_list edges;
	segment(msa, reverse_cst, node_labels, edges);
	
	std::cout << "Writing the index to disk…\n";
	make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);
	
	auto const end(chrono::high_resolution_clock::now());
	
	auto const duration(chrono::duration_cast<chrono::seconds>(end - start));
	std::cout << "Time taken: " << duration.count() << " seconds\n.";
	if (args_info.graphviz_output_given)
	{
		std::cout << "Writing the Graphviz file…\n";
		output_graphviz(node_labels, edges, args_info.graphviz_output_arg);
	}
	
	return EXIT_SUCCESS;
}
