#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/util.hpp>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <thread>         // std::thread
#include <atomic>
#include <mutex>
#include <stdio.h>
#include <unistd.h>
#include "founderblockgraph_cmdline.h"
#include "founder_block_index.hpp"

/*
   Reads MSA in fasta format as given in --input
   Finds a segment repeat-free segmentation
   Converts the segmentation into a founder block graph
   See https://arxiv.org/abs/2005.09342
   */

/* To use, install sdsl-lite: https://github.com/simongog/sdsl-lite,
   modify the Makefile, and run make */

/* Copyright (C) 2020 Veli Mäkinen under GNU General Public License v3.0 */

//#define EFG_HPP_DEBUG

namespace {

	namespace chrono = std::chrono;
	namespace fbg = founder_block_graph;

	//typedef sdsl::cst_sct3<sdsl::csa_wt<sdsl::wt_int<>>> cst_type;
	typedef sdsl::cst_sct3<> cst_type;
	typedef sdsl::csa_wt<> sa_type;
	typedef cst_type::node_type node_t;
	typedef cst_type::size_type size_type;

	typedef std::unordered_set<size_type> edge_set;
	typedef std::vector<edge_set> adjacency_list;

	inline constexpr size_type const INVALID_SIZE{std::numeric_limits <size_type>::max()};
	inline constexpr bool const VERBOSE_LOGGING{true};

#if defined(NDEBUG) && NDEBUG
	inline constexpr bool const SHOULD_STORE_ASSIGNED_NODE_LABELS{false};
#else
	inline constexpr bool const SHOULD_STORE_ASSIGNED_NODE_LABELS{true};
#endif

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


	std::string remove_gaps(std:: string const &S) {
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

	void read_input(char const *input_path, std::size_t gap_limit, bool elastic, std::vector <std::string> &msa, bool output_paths, std::vector<std::string> &identifiers)
	{
		std::string line, identifier, entry;

		// Reading input fasta
		std::fstream fs;
		fs.open(input_path, std::fstream::in);

		// Assume that the first line contains a header.
		std::getline(fs, identifier);
		if (output_paths) identifiers.push_back(identifier.substr(1));

		std::size_t expected_length(0);
		bool is_first(true);
		while (std::getline(fs, line))
		{
			if (line[0] == '>') // header
			{
				if (output_paths) identifiers.push_back(line.substr(1));
				if (is_first)
				{
					expected_length = entry.size();
					is_first = false;
				}

				if (check_sequence_length(identifier, entry, expected_length) &&
						(elastic || check_gaps(identifier, entry, gap_limit)))
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

		if (is_first)
		{
			expected_length = entry.size();
			is_first = false;
		}
		if (check_sequence_length(identifier, entry, expected_length) &&
				(elastic || check_gaps(identifier, entry, gap_limit)))
		{
			msa.push_back(entry);
		}
	} 

	void parse_input(char const *input_path, size_type &m, size_type &n, const bool output_paths, std::vector<std::string> &identifiers) {
		m = 0;
		n = 0;
		std::string line, identifier;
		std::size_t entrysize = 0;

		// Reading input fasta
		std::fstream fs;
		fs.open(input_path, std::fstream::in);

		// Assume that the first line contains a header.
		std::getline(fs, identifier);
		if (output_paths) identifiers.push_back(identifier.substr(1));

		std::size_t expected_length(0);
		bool is_first(true);
		while (std::getline(fs, line))
		{
			if (line[0] == '>') // header
			{
				if (output_paths) identifiers.push_back(line.substr(1));
				if (is_first)
				{
					expected_length = entrysize;
					n = expected_length;
					is_first = false;
				}

				if (expected_length != entrysize)
				{
					std::cerr << "MSA rows have mismatching size!" << std::endl;
					exit(1);
				}

				m += 1;
				entrysize = 0;
				identifier = line;
			}
			else
			{
				entrysize += line.size();
			}
		}
		fs.close();

		if (is_first)
		{
			expected_length = entrysize;
			n = expected_length;
		}
		if (expected_length != entrysize)
		{
			std::cerr << "MSA rows have mismatching size!" << std::endl;
			exit(1);
		}
		m += 1;
	}

	std::fstream load_rows_fs;
	std::mutex load_rows_mutex;
	std::atomic_int load_rows_startrow = 0;
	std::vector<std::string> load_rows(char const *input_path, int rows, int &output_startrow)
	{
		assert(rows > 0);
		std::vector<std::string> msa;
		std::scoped_lock lock(load_rows_mutex);
		std::string line, identifier, entry;
		output_startrow = load_rows_startrow;
		std::cerr << "Reading MSA[" << load_rows_startrow << ".." << load_rows_startrow + rows - 1 << "]..." << std::endl;

		if (!load_rows_fs.is_open())
		{
			load_rows_fs.open(input_path, std::fstream::in);
			// Assume that the first line contains a header.
			std::getline(load_rows_fs, identifier);
		}

		while (std::getline(load_rows_fs, line))
		{
			if (line[0] == '>') // header
			{
				msa.push_back(entry);
				entry.clear();
				if ((int)msa.size() >= rows) {
					load_rows_startrow += msa.size();
					return msa;
				}
			}
			else
			{
				entry += line;
			}
		}

		if (entry.size() > 0) {
			msa.push_back(entry);
		}
		load_rows_startrow += msa.size();
		return msa;
	}

	// literal C code here, TODO implement in C++ with a proper writer thread I guess
	FILE* offload_rows_fp = NULL;
	void offload_rows(char const *input_path, size_type const m, size_type const n, const std::vector<std::string> &msa, const size_type startrow)
	{
		if (offload_rows_fp == NULL) {
			offload_rows_fp = fopen((std::string(input_path) + ".transpose").c_str(), "w");
			// TODO check status
		}
		int const fd = fileno(offload_rows_fp);

		for (size_type col = 0; col < n; col++) {
			std::string strip;
			for (size_type msarow = 0; msarow < msa.size(); msarow++)
				strip += msa[msarow][col];
			pwrite(fd, strip.c_str(), msa.size(), startrow + col * m);
		}
	}

	void transpose_msa_worker(
			char const *input_path,
			size_type const m,
			size_type const n,
			const long rows
			) {
		std::string entry, identifier, line;
		std::fstream fasta_fs;
		fasta_fs.open(input_path, std::fstream::in);
		// Assume that the first line contains a header.
		std::getline(fasta_fs, identifier);

		size_type currentrow = 0;
		std::vector<std::string> msa;
		while (std::getline(fasta_fs, line))
		{
			if (line[0] == '>') // header
			{
				msa.push_back(entry);
				entry.clear();
				if ((long)msa.size() >= rows) {
					offload_rows(input_path, m, n, msa, currentrow);
					currentrow += msa.size();
				}
			}
			else
			{
				entry += line;
			}
		}

		if (entry.size() > 0) {
			msa.push_back(entry);
		}
		offload_rows(input_path, m, n, msa, currentrow);
		currentrow += msa.size();
		fasta_fs.close();
	}

	bool load_cst(char const *input_path, std::vector<std::string> const &MSA, cst_type &cst, size_t gap_limit) {

		std::string const index_suffix(".cst");
		std::string const plain_suffix(".plain");
		std::string index_file = std::string(input_path) + plain_suffix + std::to_string(gap_limit) + index_suffix;

		// Constructing compressed suffix tree for C
		//if (!sdsl::load_from_file(cst, index_file))
		{
			std::cerr << "No index "<<index_file<< " located. Building index now." << std::endl;

			// Output concatenated inputs to disk. Remove gap symbols and add separators for indexing.
			{
				std::fstream fs;
				fs.open(std::string(input_path) + plain_suffix, std::fstream::out);
				for (auto const &seq : MSA)
				{
					for (auto const c : seq)
					{
						if ('-' != c)
							fs << c;
					}
					fs << '#';
				}
				fs.close();
			}
			/*{
				size_type size = 0;
				for (auto const &seq : MSA) {
					size += remove_gaps(seq).size() + 1;
				}
				sdsl::int_vector<> todisk(size, 0, 4);
				size_type i = 0;
				for (auto const &seq : MSA) {
					for (auto const c : seq) {
						if (c == '-')
							continue;

						if (c == 'A') {
							todisk[i++] = 1;
						} else if (c == 'C') {
							todisk[i++] = 2;
						} else if (c == 'G') {
							todisk[i++] = 3;
						} else if (c == 'T') {
							todisk[i++] = 4;
						} else {
							todisk[i++] = 5;
						}
					}
					// separator character
					todisk[i++] = 6;
				}
				//sdsl::util::bit_compress(todisk);
				std::fstream fs;
				//fs.open(std::string(input_path) + plain_suffix, std::fstream::out);
				//todisk.serialize(fs);
				//fs.close();
				sdsl::store_to_file(todisk, std::string(input_path) + plain_suffix);
			}*/

			std::ifstream in(std::string(input_path)+plain_suffix);
			if (!in) {
				std::cerr << "ERROR: File " << input_path << ".plain" << " does not exist. Exit." << std::endl;
				return false;
			}
			//sdsl::construct(cst, std::string(input_path)+plain_suffix, 0); // generate index
			sdsl::construct(cst, std::string(input_path)+plain_suffix, 1); // generate index
			//sdsl::csXprintf(std::cout, "%3I%3S %3s %3P %3p %3L %3B  %T", cst, '$');
			sdsl::store_to_file(cst, index_file); // save it
		}/* else {
			std::cerr << "Index "<<index_file<< " located. Did you use option --heuristic-subset? If so, delete the index." << std::endl;
		}*/

		return true;
	}
	bool compute_cst_im(std::vector<std::string> const &MSA, cst_type &cst) {
		std::string msaconcat;

		// Constructing compressed suffix tree for C
		//if (!sdsl::load_from_file(cst, index_file))
		{
			// compute concatenated inputs. Remove gap symbols and add separators for indexing.
			{
				for (auto const &seq : MSA)
				{
					for (auto const c : seq)
					{
						if ('-' != c)
							msaconcat += c;
					}
					msaconcat += '#';
				}
			}

			//sdsl::construct(cst, std::string(input_path)+plain_suffix, 0); // generate index
			sdsl::construct_im(cst, msaconcat, 1); // generate index
			//sdsl::csXprintf(std::cout, "%3I%3S %3s %3P %3p %3L %3B  %T", cst, '$');
			//sdsl::store_to_file(cst, index_file); // save it
		}/* else {
			std::cerr << "Index "<<index_file<< " located. Did you use option --heuristic-subset? If so, delete the index." << std::endl;
		}*/

		return true;
	}
	bool load_sa(char const *input_path, std::vector<std::string> const &MSA, sa_type &sa) {

		std::string const index_suffix(".sa");
		std::string const plain_suffix(".sa.plain");
		std::string index_file = std::string(input_path) + plain_suffix + index_suffix;

		// Constructing compressed suffix tree for C
		if (!sdsl::load_from_file(sa, index_file))
		{
			std::cerr << "No index "<<index_file<< " located. Building index now." << std::endl;

			// Output concatenated inputs to disk. Remove gap symbols and add separators for indexing.
			{
				std::fstream fs;
				fs.open(std::string(input_path) + plain_suffix, std::fstream::out);
				for (auto const &seq : MSA)
				{
					for (unsigned long i = 0; i < seq.size(); i++)
					{
						if ('-' != seq[i])
							fs << seq[i];
					}
					fs << '#';
				}
				fs.close();
			}

			std::ifstream in(std::string(input_path)+plain_suffix);
			if (!in) {
				std::cerr << "ERROR: File " << input_path << ".plain" << " does not exist. Exit." << std::endl;
				return false;
			}
			sdsl::construct(sa, std::string(input_path)+plain_suffix, 1); // generate index
			sdsl::store_to_file(sa, index_file); // save it
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

	bool interval::operator<(interval const &other) const
	{
		if (sp < other.sp) return 1;
		else if (sp == other.sp && ep > other.ep) return 1; // *this contains other but is not equivalent with it.
		else return 0;
	}

	int segment(
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::vector <std::string> &out_labels,
			adjacency_list &out_edges
		   ) {

		size_type n = MSA[0].size();

		/* v[j] is such that MSA[0..m-1][v[j]..j] is a repeat-free segment */
		/* This is a preprocessing part of the main algorithm */

		std::vector<size_type> v(n, 0);

		size_type m = MSA.size();
		std::vector<size_type> sp(m, 0);
		std::vector<size_type> ep(m, cst.csa.size()-1);
		std::vector<size_type> lcp(m, 0);
		std::vector<interval> pairs;
		/* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
		/* lcp[i] = j-jp+1 minus the number of gap symbols in that range*/
		/* pairs contains (sp[i],ep[i]) and is used for sorting the intervals */

		// TODO: replace the computation of v[j]'s with call to find_valid_blocks()     
		size_type jp = n; // maintains the left-boundary 
		size_type sum;
		for (size_type j=n-1; j+1>0; j--) {
			v[j] = j+1; // no valid range found
			if (j<n-1) {  
				/* contracting MSA[i][j] from right for all i */
				for (size_type i=0; i<m; i++) {
					/* mapping BWT range to suffix tree node,
					   taking parent, mapping back if parent is not too far*/
					// Keeping old range on gap symbols
					if (MSA[i][j+1]!='-') {    
						node_t ll = cst.select_leaf(sp[i]+1); 
						node_t rl = cst.select_leaf(ep[i]+1); 
						node_t w = cst.lca(ll,rl);
						node_t u = cst.parent(w);
						lcp[i]--;              
						if (cst.depth(u)==lcp[i]) {
							// contracting first symbol from edge (u,w)
							sp[i] = cst.lb(u);
							ep[i] = cst.rb(u);
						}
					} 
				}
			}   
			while (1) {
				/* Now checking if the union of intervals is of size m */
				/* If that is the case, we have found v(j), otherwise continue left-extensions */
				/* Intervals are distinct or nested, proper overlaps are impossible */
				/* We can thus sort primarily by smallest left-boundary and secondary by largest right-boundary*/     
				/* Then intervals get nested, and we can just sum the sizes of the out-most intervals*/     
				pairs.clear();
				for (size_type i=0; i<m; i++)            
					pairs.push_back(interval(sp[i],ep[i]));
				std::sort(pairs.begin(),pairs.end());
				sum = 0;
				size_type spprev = pairs[0].sp;
				size_type epprev = pairs[0].ep;
				for (size_type i=1; i<m; i++)
					if (pairs[i].sp>epprev) {
						sum += epprev-spprev+1;
						spprev= pairs[i].sp; 
						epprev = pairs[i].ep;          
					}   
				sum += epprev-spprev+1;  
				if (sum == m) { // no non-aligned matches
					v[j]=jp;
					break;
				}    
				if (jp==0) // cannot go more to the left
					break;        
				jp--; 
				/* Now looking at MSA[0..m-1][jp..j] */
				for (size_type i=0; i<m; i++) {
					// Keeping old range on gap symbols
					assert(jp < MSA[i].size());
					if (MSA[i][jp]!='-') {
						sdsl::backward_search(cst.csa,sp[i],ep[i],MSA[i][jp],sp[i],ep[i]);
						lcp[i]++;
					}
				}
			} 
		}  

		/* Main algorithm begins */
		/* s[j] is the score of the optimal valid segmentation of MSA[0..m-1][0..j] */
		/* s[n-1]=min_{S is valid segmentation} max_{[a..b] \in S} b-a+1 */ 
		std::vector<size_type> s(n, n);

		/* prev[j] is the pointer to the beginning of the segment ending at j in optimal segmentation */
		std::vector<size_type> prev(n, n);

		// Computation of s[j]'s and prev[j]'s
		for (size_type j=0; j<n; j++) {
			s[j] = j+2; // no valid range
			prev[j] = j+1; // no valid range
			if (v[j]>j) 
				continue; // no valid range
			jp=v[j];  
			while (1) {
				if (jp!=0 && s[jp-1]==jp+1)  {// no proper segmentation ending at jp-1
					jp--;
					continue; 
				}
				if (s[j]>std::max(jp==0?0:s[jp-1],j-jp+1)) {
					s[j] = std::max(jp==0?0:s[jp-1],j-jp+1);
					prev[j] = jp;       
				}
				if (s[j]==j-jp+1) 
					break;
				if (jp==0)
					break;
				jp--;
			}
		}

		// outputing optimal score 
		std::cerr << "Optimal score: " << s[n-1] << std::endl;

		if (s[n-1]==n+1) // no proper segmentation exists
		{
			std::cerr << "No proper segmentation exists.\n";
			return EXIT_FAILURE;
		}

		std::list<size_type> boundariestemp;
		size_type j = n-1;
		boundariestemp.push_front(j);
		while (prev[j]!=0) {
			boundariestemp.push_front(prev[j]-1);
			j = prev[j]-1;         
		}
		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);
		std::cerr << "Number of segments: " << boundaries.size() << std::endl;

		/* Convert segmentation into founder block graph */
		std::unordered_map<std::string, size_type> str2id;
		size_type nodecount = 0;
		size_type previndex = 0;

		typedef std::vector<size_type> block_vector;
		typedef std::vector<block_vector> block_matrix;
		block_matrix blocks(boundaries.size());
		std::string ellv, ellw;
		for (size_type j=0; j<boundaries.size(); j++) {
			for (size_type i=0; i<m; i++) {
				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));   
				if (!str2id.count(ellv)) {  
					blocks[j].push_back(nodecount);
					str2id[ellv] = nodecount++;
				}
			}
			previndex = boundaries[j]+1;
		}

		std::vector<std::string> labels(nodecount);
		for (const auto& pair : str2id) {
			labels[pair.second] = pair.first;
		}
		size_type totallength = 0;
		for (size_type i=0; i<nodecount; i++)
			totallength += labels[i].size();

		std::cerr << "#nodes=" << nodecount << std::endl;
		std::cerr << "total length of node labels=" << totallength << std::endl;

		size_type nfounders=0;
		for (size_type j=0; j<boundaries.size(); j++)
			if (blocks[j].size()>nfounders)
				nfounders = blocks[j].size();
		std::cerr << "#founders=" << nfounders << std::endl;

		adjacency_list edges(nodecount);
		previndex = 0;
		for (size_type k=0; k<boundaries.size()-1; k++)
		{
			for (size_type i=0; i<m; i++)
			{

				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[k]-previndex+1));
				ellw = remove_gaps(MSA[i].substr(boundaries[k]+1,boundaries[k+1]-boundaries[k]));
				auto const &src_node_label(ellv);
				auto const &dst_node_label(ellw);
				auto const src_node_idx_it(str2id.find(src_node_label));
				auto const dst_node_idx_it(str2id.find(dst_node_label));
				assert(src_node_idx_it != str2id.end());
				assert(dst_node_idx_it != str2id.end());
				auto const src_node_idx(src_node_idx_it->second);
				auto const dst_node_idx(dst_node_idx_it->second);
				edges[src_node_idx].insert(dst_node_idx);
			}

			previndex = boundaries[k]+1;
		}
		size_type edgecount = 0;
		for (size_type i=0; i<nodecount; i++)
			edgecount += edges[i].size();
		std::cerr << "#edges=" << edgecount << std::endl;

		using std::swap;
		swap(out_labels, labels);
		swap(out_edges, edges);

		return EXIT_SUCCESS;
	}


	int segment2elasticValid(
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::vector <std::string> &out_labels,
			adjacency_list &out_edges
			) {

		size_type n = MSA[0].size();

		size_type m = MSA.size();

		/* v[j] is such that MSA[0..m-1][v[j]..j] is a repeat-free segment */
		/* This is a preprocessing part of the main algorithm */

		std::vector<size_type> v(n, 0);

		std::vector<size_type> sp(m, 0);
		std::vector<size_type> ep(m, cst.csa.size()-1);
		std::vector<size_type> lcp(m, 0);
		std::vector<interval> pairs;
		/* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
		/* lcp[i] = j-jp+1 minus the number of gap symbols in that range*/
		/* pairs contains (sp[i],ep[i]) and is used for sorting the intervals */
		size_type jp = n; // maintains the left-boundary 
		size_type sum;
		for (size_type j=n-1; j+1>0; j--) {
			v[j] = j+1; // no valid range found
			if (j<n-1) {  
				/* contracting MSA[i][j] from right for all i */
				for (size_type i=0; i<m; i++) {
					/* mapping BWT range to suffix tree node,
					   taking parent, mapping back if parent is not too far*/
					// Keeping old range on gap symbols
					if (MSA[i][j+1]!='-') {    
						node_t ll = cst.select_leaf(sp[i]+1); 
						node_t rl = cst.select_leaf(ep[i]+1); 
						node_t w = cst.lca(ll,rl);
						node_t u = cst.parent(w);
						lcp[i]--;              
						if (cst.depth(u)==lcp[i]) {
							// contracting first symbol from edge (u,w)
							sp[i] = cst.lb(u);
							ep[i] = cst.rb(u);
						}
					} 
				}
			}   
			while (1) {
				/* Now checking if the union of intervals is of size m */
				/* If that is the case, we have found v(j), otherwise continue left-extensions */
				/* Intervals are distinct or nested, proper overlaps are impossible */
				/* We can thus sort primarily by smallest left-boundary and secondary by largest right-boundary*/     
				/* Then intervals get nested, and we can just sum the sizes of the out-most intervals*/     
				pairs.clear();
				for (size_type i=0; i<m; i++)            
					pairs.push_back(interval(sp[i],ep[i]));
				std::sort(pairs.begin(),pairs.end());
				sum = 0;
				size_type spprev = pairs[0].sp;
				size_type epprev = pairs[0].ep;
				for (size_type i=1; i<m; i++)
					if (pairs[i].sp>epprev) {
						sum += epprev-spprev+1;
						spprev= pairs[i].sp; 
						epprev = pairs[i].ep;          
					}   
				sum += epprev-spprev+1;  
				if (sum == m) { // no non-aligned matches
					v[j]=jp;
					break;
				}    
				if (jp==0) // cannot go more to the left
					break;        
				jp--; 
				/* Now looking at MSA[0..m-1][jp..j] */
				for (size_type i=0; i<m; i++) {
					// Keeping old range on gap symbols
					assert(jp < MSA[i].size());
					if (MSA[i][jp]!='-') {
						sdsl::backward_search(cst.csa,sp[i],ep[i],MSA[i][jp],sp[i],ep[i]);
						lcp[i]++;
					}
				}
			} 
		}  
		/* Main algorithm begins */
		/* s[j] is the score of the valid elastic segmentation of MSA[0..m-1][0..j] */
		/* s[n-1]=min_{S is valid elastic segmentation} max_{[a..b] \in S} b-a+1 */
		/* Heuristic applied so that s[n-1] is valid but not necessarily optimal */
		std::vector<size_type> s(n, n+1);
		/* prev[j] is the pointer to the end of previous segment in optimal segmentation */
		std::vector<size_type> prev(n, n+1);
		for (size_type j=1; j<n; j++) {
			jp = v[j]; // start of valid segment
			if (jp>j) 
				continue;
			else if (jp==0) {
				s[j] = j+1;
				prev[j] = 0;  
			}
			else if (std::max(s[jp-1],j-jp+1)<std::max(s[j-1],j-prev[j-1]+1)) {
				s[j] = std::max(s[jp-1],j-jp+1);
				prev[j] = jp;
			}
			else {
				s[j] = std::max(s[j-1],j-prev[j-1]+1);
				prev[j] = prev[j-1];  
			}
		}  
		// outputing optimal score 
		std::cerr << "Optimal score: " << s[n-1] << std::endl;

		if (s[n-1]==n+1) // no proper segmentation exists
		{
			std::cerr << "No valid segmentation found!\n";
			return EXIT_FAILURE;
		}

		std::list<size_type> boundariestemp;
		size_type j = n-1;
		boundariestemp.push_front(j);
		while (prev[j]!=0) {
			boundariestemp.push_front(prev[j]-1);
			j = prev[j]-1;         
		}
		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);
		std::cerr << "Number of segments: " << boundaries.size() << std::endl;

		/* Convert segmentation into founder block graph */
		std::unordered_map<std::string, size_type> str2id;
		size_type nodecount = 0; 
		size_type previndex = 0;

		typedef std::vector<size_type> block_vector;
		typedef std::vector<block_vector> block_matrix;
		block_matrix blocks(boundaries.size());
		std::string ellv, ellw;
		for (size_type j=0; j<boundaries.size(); j++) {
			for (size_type i=0; i<m; i++) {
				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));
				if (!str2id.count(ellv)) {       
					blocks[j].push_back(nodecount);
					str2id[ellv] = nodecount++;
				}
			}
			previndex = boundaries[j]+1;
		}

		std::vector<std::string> labels(nodecount);
		for (const auto& pair : str2id) {
			labels[pair.second] = pair.first;
		}
		size_type totallength = 0;
		for (size_type i=0; i<nodecount; i++)
			totallength += labels[i].size();

		std::cerr << "#nodes=" << nodecount << std::endl;
		std::cerr << "total length of node labels=" << totallength << std::endl;

		size_type nfounders=0;
		for (size_type j=0; j<boundaries.size(); j++)
			if (blocks[j].size()>nfounders)
				nfounders = blocks[j].size();
		std::cerr << "#founders=" << nfounders << std::endl;

		adjacency_list edges(nodecount);
		previndex = 0;
		for (size_type k=0; k<boundaries.size()-1; k++)
		{
			for (size_type i=0; i<m; i++)
			{
				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[k]-previndex+1));
				ellw = remove_gaps(MSA[i].substr(boundaries[k]+1,boundaries[k+1]-boundaries[k]));
				auto const &src_node_label(ellv);
				auto const &dst_node_label(ellw);
				auto const src_node_idx_it(str2id.find(src_node_label));
				auto const dst_node_idx_it(str2id.find(dst_node_label));
				assert(src_node_idx_it != str2id.end());
				assert(dst_node_idx_it != str2id.end());
				auto const src_node_idx(src_node_idx_it->second);
				auto const dst_node_idx(dst_node_idx_it->second);
				edges[src_node_idx].insert(dst_node_idx);
			}

			previndex = boundaries[k]+1;
		}
		size_type edgecount = 0;
		for (size_type i=0; i<nodecount; i++)
			edgecount += edges[i].size();
		std::cerr << "#edges=" << edgecount << std::endl;

		using std::swap;
		swap(out_labels, labels);
		swap(out_edges, edges);
		return EXIT_SUCCESS;
	}

	void make_efg(
			const std::vector<size_type> &boundaries,
			const std::vector<std::string> &MSA,
			std::vector <std::string> &out_labels,
			std::vector <size_type> &out_blocks,
			adjacency_list &out_edges,
			bool output_paths,
			std::vector<std::vector<size_type>> &out_paths
	) {
		out_labels.clear();
		out_blocks.clear();
		out_edges.clear();
		out_paths.clear();

		/* Convert arbitrary segmentation into EFG */
		// TODO: check, O(log(n)) complexity of operations?
		size_type const m = MSA.size();

		//std::unordered_map<std::string, std::pair<size_type,size_type>> str2id; // map label -> node_id, block?
		std::vector<std::unordered_map<std::string, std::pair<size_type,size_type>>> str2ids (boundaries.size()); // map label -> node_id, block?
		size_type nodecount = 0;
		size_type previndex = 0;

		typedef std::vector<size_type> block_vector;
		typedef std::vector<block_vector> block_matrix;
		block_matrix blocks(boundaries.size());
		std::string ellv, ellw;
		std::vector<std::vector<size_type>> paths (m, std::vector<size_type>());
		for (size_type j=0; j<boundaries.size(); j++) {
			for (size_type i=0; i<m; i++) {
				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));
				if (ellv.length() == 0)
					continue;
				if (!str2ids[j].count(ellv)) {
					blocks[j].push_back(nodecount);
					std::pair p(nodecount++, j);
					str2ids[j][ellv] = p;
					//std::cerr << "p.first = " << p.first << " and p.second = " << p.second << "\n";
				}
				if (output_paths) paths[i].push_back(str2ids[j][ellv].first);
			}
			previndex = boundaries[j]+1;
		}
		/*std::cerr << "My paths are: " << std::endl;
		  for (size_type i = 0; i < m; i++) {
		  std::cerr << i << ": ";
		  for (auto s : paths[i]) {
		  std::cerr << s+1 << " ";
		  }
		  std::cerr << std::endl;
		  }*/

		std::vector<std::string> labels(nodecount);
		std::vector<size_type> bblocks(nodecount);
		for (auto &str2id : str2ids) {
			for (const auto& pair : str2id) {
				labels[pair.second.first] = pair.first;
				bblocks[pair.second.first] = pair.second.second;
			}
		}
		size_type totallength = 0;
		for (size_type i=0; i<nodecount; i++)
			totallength += labels[i].size();

		std::cerr << "#nodes=" << nodecount << std::endl;
		std::cerr << "total length of node labels=" << totallength << std::endl;

		/*size_type nfounders=0;
		  for (size_type j=0; j<boundaries.size(); j++)
		  if (blocks[j].size()>nfounders)
		  nfounders = blocks[j].size();
		  std::cerr << "#founders=" << nfounders << std::endl;*/

		adjacency_list edges(nodecount);
		previndex = 0;
		for (size_type k=0; k<boundaries.size()-1; k++)
		{
			for (size_type i=0; i<m; i++)
			{
				ellv = remove_gaps(MSA[i].substr(previndex,boundaries[k]-previndex+1));
				ellw = remove_gaps(MSA[i].substr(boundaries[k]+1,boundaries[k+1]-boundaries[k]));
				if (ellv.length() == 0 || ellw.length() == 0)
					continue;
				auto const &src_node_label(ellv);
				auto const &dst_node_label(ellw);
				auto const src_node_idx_it(str2ids[k].find(src_node_label));
				auto const dst_node_idx_it(str2ids[k+1].find(dst_node_label));
				assert(src_node_idx_it != str2ids[k].end());
				assert(dst_node_idx_it != str2ids[k+1].end());
				auto const src_node_idx(src_node_idx_it->second);
				auto const dst_node_idx(dst_node_idx_it->second);
				edges[src_node_idx.first].insert(dst_node_idx.first);
			}

			previndex = boundaries[k]+1;
		}
		size_type edgecount = 0;
		for (size_type i=0; i<nodecount; i++)
			edgecount += edges[i].size();
		std::cerr << "#edges=" << edgecount << std::endl;

		using std::swap;
		swap(out_labels, labels);
		swap(out_blocks, bblocks);
		swap(out_edges, edges);
		swap(out_paths, paths);
		std::cerr << std::flush;
	}

	void make_efg_external(
			char const *input_path,
			size_type m,
			size_type n,
			const std::vector<size_type> &boundaries,
			const std::vector<std::string> &MSA,
			std::vector <std::string> &out_labels,
			std::vector <size_type> &out_blocks,
			adjacency_list &out_edges,
			bool output_paths,
			std::vector<std::vector<size_type>> &out_paths
	) {
		std::fstream fs;
		fs.open(std::string(input_path) + ".transpose", std::fstream::in);

		out_labels.clear();
		out_blocks.clear();
		out_edges.clear();
		out_paths.clear();

		/* Convert arbitrary segmentation into EFG */
		// TODO: check, O(log(n)) complexity of operations?

		std::vector<std::unordered_map<std::string, std::pair<size_type,size_type>>> str2ids (boundaries.size()); // map label -> node_id, block?
		size_type nodecount = 0;
		size_type previndex = 0;

		typedef std::vector<size_type> block_vector;
		typedef std::vector<block_vector> block_matrix;
		block_matrix blocks(boundaries.size());
		std::string ellv, ellw;
		std::vector<std::vector<size_type>> paths(m, std::vector<size_type>());
		for (size_type j=0; j<boundaries.size(); j++) {
			// read (boundaries[j] - previndex + 1) columns
			std::vector<std::string> block_transpose;
			if (boundaries[j]+1 <= previndex) {
				std::cerr << "DEBUG: boundaries[j] is " << boundaries[j] << " and previndex is " << previndex << std::endl;
			}
			assert(boundaries[j]+1 > previndex);
			for (size_type i = 0; i < boundaries[j] - previndex + 1; i++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose.push_back(std::string(s.data()));
			}
			for (size_type i=0; i<m; i++) {
				ellv.clear();
				for (size_type j = 0; j < block_transpose.size(); j++)
					ellv += block_transpose[j][i];
				ellv = remove_gaps(ellv);
				if (ellv.length() == 0)
					continue;
				if (!str2ids[j].count(ellv)) {
					blocks[j].push_back(nodecount);
					std::pair p(nodecount++, j);
					str2ids[j][ellv] = p;
				}
				if (output_paths) paths[i].push_back(str2ids[j][ellv].first);
			}
			previndex = boundaries[j]+1;
		}

		std::vector<std::string> labels(nodecount);
		std::vector<size_type> bblocks(nodecount);
		for (auto &str2id : str2ids) {
			for (const auto& pair : str2id) {
				labels[pair.second.first] = pair.first;
				bblocks[pair.second.first] = pair.second.second;
			}
		}
		size_type totallength = 0;
		for (size_type i=0; i<nodecount; i++)
			totallength += labels[i].size();

		std::cerr << "#nodes=" << nodecount << std::endl;
		std::cerr << "total length of node labels=" << totallength << std::endl;

		adjacency_list edges(nodecount);
		previndex = 0;
		for (size_type k=0; k<boundaries.size()-1; k++)
		{
			std::vector<std::string> block_transpose, block_transpose2;
			fs.clear();
			fs.seekg(previndex * m);
			for (size_type ii = 0; ii < boundaries[k]-previndex+1; ii++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose.push_back(std::string(s.data()));
			}
			fs.clear();
			fs.seekg((boundaries[k]+1) * m);
			for (size_type ii = 0; ii < boundaries[k+1]-boundaries[k]; ii++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose2.push_back(std::string(s.data()));
			}
			for (size_type i=0; i<m; i++)
			{
				ellv.clear();
				for (size_type j = 0; j < block_transpose.size(); j++)
					ellv += block_transpose[j][i];
				ellv = remove_gaps(ellv);

				ellw.clear();
				for (size_type j = 0; j < block_transpose2.size(); j++)
					ellw += block_transpose2[j][i];
				ellw = remove_gaps(ellw);

				if (ellv.length() == 0 || ellw.length() == 0)
					continue;
				auto const &src_node_label(ellv);
				auto const &dst_node_label(ellw);
				auto const src_node_idx_it(str2ids[k].find(src_node_label));
				auto const dst_node_idx_it(str2ids[k+1].find(dst_node_label));
				assert(src_node_idx_it != str2ids[k].end());
				assert(dst_node_idx_it != str2ids[k+1].end());
				auto const src_node_idx(src_node_idx_it->second);
				auto const dst_node_idx(dst_node_idx_it->second);
				edges[src_node_idx.first].insert(dst_node_idx.first);
			}

			previndex = boundaries[k]+1;
		}
		size_type edgecount = 0;
		for (size_type i=0; i<nodecount; i++)
			edgecount += edges[i].size();
		std::cerr << "#edges=" << edgecount << std::endl;

		using std::swap;
		swap(out_labels, labels);
		swap(out_blocks, bblocks);
		swap(out_edges, edges);
		swap(out_paths, paths);
		std::cerr << std::flush;
		fs.close();
	}

	void output_efg(
			const std::vector<size_type> &boundaries,
			const std::vector<std::string> &MSA,
			bool output_paths,
			std::vector<std::string> const &identifiers,
			char const *dst_path_c
	) {
		size_type const n = MSA[0].size();
		size_type const m = MSA.size();
		// create the output file
		std::string const dst_path(dst_path_c);
		std::fstream dst_os;
		dst_os.exceptions(std::fstream::failbit);
		dst_os.open(dst_path, std::ios_base::out);

		// MSA info
		dst_os << "M\t" << m << "\t" << n << std::endl;

		// Segmentation info
		dst_os << "X\t1";
		for (size_type i = 0; (int)i < (int)boundaries.size() - 1; i++)
			dst_os << "\t" << boundaries[i] + 2;
		dst_os << std::endl;

		// Partition of nodes into blocks
		dst_os << "B\t";
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_set<std::string> labels;
			for (size_type i=0; i<m; i++) {
				std::string const label = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));
				if (label != "")
					labels.insert(label);
			}
			dst_os << ((j == 0) ? "" : "\t") << labels.size();
		}
		dst_os << std::endl;

		// Output nodes and edges
		// data structures to map node labels of a block to their ID/index
		std::unordered_map<std::string, size_type> str2id_previous_block;
		std::unordered_map<size_type,   size_type> row2id_previous_block;
		std::unordered_map<std::string, size_type> str2id_current_block;
		size_type nodecount = 0;
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_map<size_type, size_type> row2id_current_block;
			std::set<std::pair<size_type, size_type>> edges_to_previous_block;

			for (size_type i=0; i<m; i++) {
				// node
				std::string label = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));
				if (label == "") // empty label
					continue;

				size_type nodeindex;
				if (str2id_current_block.find(label) == str2id_current_block.end()) {
					nodeindex = nodecount++;
					dst_os << "S\t" << nodeindex << "\t" << label << std::endl;
				} else {
					nodeindex = str2id_current_block[label];
				}
				str2id_current_block[label] = nodeindex;
				row2id_current_block[i] = nodeindex;

				// compute edge
				if (row2id_previous_block.find(i) != row2id_previous_block.end())
					edges_to_previous_block.insert(std::pair<size_type, size_type>(row2id_previous_block[i], nodeindex));
			}
			// print edges
			for (auto &p : edges_to_previous_block) {
				dst_os << "L\t" << p.first << "\t+\t" << p.second << "\t+\t0M" << std::endl;
			}

			swap(row2id_previous_block, row2id_current_block);
			swap(str2id_previous_block, str2id_current_block);
			str2id_current_block.clear();
		}

		// we are done, if we do not need to output the paths
		if (!output_paths)
			return;

		// TODO: int_vector of integers from 0 to H
		std::vector<std::vector<size_type>> paths (m, std::vector<size_type>());
		nodecount = 0;
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_map<std::string, size_type> str2id_current_block;
			std::unordered_map<size_type,   size_type> row2id_current_block;
			for (size_type i=0; i<m; i++) {
				std::string const label = remove_gaps(MSA[i].substr(previndex,boundaries[j]-previndex+1));
				if (label == "") // empty label
					continue;

				size_type nodeindex;
				if (str2id_current_block.find(label) == str2id_current_block.end()) {
					nodeindex = nodecount++;
				} else {
					nodeindex = str2id_current_block[label];
				}
				str2id_current_block[label] = nodeindex;
				row2id_current_block[i] = nodeindex;
			}
			for (auto &p : row2id_current_block) {
				paths[p.first].push_back(p.second);
			}
		}

		// paths
		assert(identifiers.size() == paths.size());
		for (size_type i = 0; i < paths.size(); i++) {
			dst_os << "P\t" << identifiers[i] << "\t";
			for (size_type j = 0; j < paths[i].size() - 1; j++) {
				dst_os << paths[i][j] << "+,";
			}
			dst_os << paths[i][paths[i].size()-1] << "+";
			dst_os << "\t*" << std::endl;
		}
	}

	void output_efg_external(
			const std::vector<size_type> &boundaries,
			char const *input_path,
			size_type const m,
			size_type const n,
			bool output_paths,
			std::vector<std::string> const &identifiers,
			char const *dst_path_c
	) {
		std::fstream fs;
		fs.open(std::string(input_path) + ".transpose", std::fstream::in);

		// create the output file
		std::string const dst_path(dst_path_c);
		std::fstream dst_os;
		dst_os.exceptions(std::fstream::failbit);
		dst_os.open(dst_path, std::ios_base::out);

		// MSA info
		dst_os << "M\t" << m << "\t" << n << std::endl;

		// Segmentation info
		dst_os << "X\t1";
		for (size_type i = 0; (int)i < (int)boundaries.size() - 1; i++)
			dst_os << "\t" << boundaries[i] + 2;
		dst_os << std::endl;

		// Partition of nodes into blocks
		dst_os << "B\t";
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_set<std::string> labels;
			std::vector<std::string> block_transpose;
			for (size_type i = 0; i < boundaries[j] - previndex + 1; i++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose.push_back(std::string(s.data()));
			}

			for (size_type i=0; i<m; i++) {
				std::string label;
				for (size_type j = 0; j < block_transpose.size(); j++)
					label += block_transpose[j][i];
				label = remove_gaps(label);
				if (label != "")
					labels.insert(label);
			}
			dst_os << ((j == 0) ? "" : "\t") << labels.size();
		}
		dst_os << std::endl;

		// Output nodes and edges
		// data structures to map node labels of a block to their ID/index
		std::unordered_map<std::string, size_type> str2id_previous_block;
		std::unordered_map<size_type,   size_type> row2id_previous_block;
		std::unordered_map<std::string, size_type> str2id_current_block;
		size_type nodecount = 0;
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_map<size_type, size_type> row2id_current_block;
			std::set<std::pair<size_type, size_type>> edges_to_previous_block;

			std::vector<std::string> block_transpose;
			fs.clear();
			fs.seekg(previndex * m);
			for (size_type ii = 0; ii < boundaries[j]-previndex+1; ii++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose.push_back(std::string(s.data()));
			}

			for (size_type i=0; i<m; i++) {
				// node
				std::string label;
				for (size_type j = 0; j < block_transpose.size(); j++)
					label += block_transpose[j][i];
				label = remove_gaps(label);
				if (label == "") // empty label
					continue;

				size_type nodeindex;
				if (str2id_current_block.find(label) == str2id_current_block.end()) {
					nodeindex = nodecount++;
					dst_os << "S\t" << nodeindex << "\t" << label << std::endl;
				} else {
					nodeindex = str2id_current_block[label];
				}
				str2id_current_block[label] = nodeindex;
				row2id_current_block[i] = nodeindex;

				// compute edge
				if (row2id_previous_block.find(i) != row2id_previous_block.end())
					edges_to_previous_block.insert(std::pair<size_type, size_type>(row2id_previous_block[i], nodeindex));
			}
			// print edges
			for (auto &p : edges_to_previous_block) {
				dst_os << "L\t" << p.first << "\t+\t" << p.second << "\t+\t0M" << std::endl;
			}

			swap(row2id_previous_block, row2id_current_block);
			swap(str2id_previous_block, str2id_current_block);
			str2id_current_block.clear();
		}

		// we are done, if we do not need to output the paths
		if (!output_paths)
			return;

		// TODO: int_vector of integers from 0 to H
		std::vector<std::vector<size_type>> paths (m, std::vector<size_type>());
		nodecount = 0;
		for (size_t j = 0, previndex = 0; (int)j < (int)boundaries.size(); previndex = boundaries[j]+1, j++) {
			std::unordered_map<std::string, size_type> str2id_current_block;
			std::unordered_map<size_type,   size_type> row2id_current_block;

			std::vector<std::string> block_transpose;
			fs.clear();
			fs.seekg(previndex * m);
			assert(boundaries[j] > previndex - 1);
			for (size_type ii = 0; ii < boundaries[j]-previndex+1; ii++) {
				std::vector<char> s(m+1);
				fs.get(&s[0], m+1);
				if (s[0] != 0)
					block_transpose.push_back(std::string(s.data()));
			}

			for (size_type i=0; i<m; i++) {
				std::string label;
				for (size_type j = 0; j < block_transpose.size(); j++)
					label += block_transpose[j][i];
				label = remove_gaps(label);
				if (label == "") // empty label
					continue;

				size_type nodeindex;
				if (str2id_current_block.find(label) == str2id_current_block.end()) {
					nodeindex = nodecount++;
				} else {
					nodeindex = str2id_current_block[label];
				}
				str2id_current_block[label] = nodeindex;
				row2id_current_block[i] = nodeindex;
			}
			for (auto &p : row2id_current_block) {
				paths[p.first].push_back(p.second);
			}
		}

		// paths
		assert(identifiers.size() == paths.size());
		for (size_type i = 0; i < paths.size(); i++) {
			dst_os << "P\t" << identifiers[i] << "\t";
			for (size_type j = 0; j < paths[i].size() - 1; j++) {
				dst_os << paths[i][j] << "+,";
			}
			dst_os << paths[i][paths[i].size()-1] << "+";
			dst_os << "\t*" << std::endl;
		}

		fs.close();
	}

	/* 
	 * compute_f_range accepts in input the MSA, its compressed suffix tree,
	 * the relative sdsl data structures on its rows and a column interval
	 * [startx..endx]. It computes the minimal right extension of x for x
	 * in [startx..endx] and stores the values in f[x], such that
	 * f[x] in minimum index >= x such that MSA[0..m-1][x..f[x]] is
	 * semi-repeat-free.
	 *
	 * This implementation is O(l * m * log m), where l = endx - startx.
	 */
	void compute_f_range(
			size_type const m,
			size_type const n,
			size_type const startx,
			size_type const endx,
			std::vector<size_type> &f, // store result here
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::string const &ignorechars,
			std::vector<sdsl::rank_support_v5<>> const &rowsignore_rs,
			std::vector<sdsl::select_support_mcl<>> const &rowsignore_ss,
			sdsl::rank_support_v5<> const &concatenated_rows_rs,
			std::vector<sdsl::rank_support_v5<>> const &indexedrows_rs,
			std::vector<sdsl::select_support_mcl<>> const &indexedrows_ss,
			const bool disable_efg_tricks
			) {
		// find leaves corresponding to suffixes starting at col x+1
		std::vector<node_t> leaves(m, cst.root());
		std::unordered_map<size_type, size_type>  leavesmap; // leaf index -> MSA row
		for (size_type indexpos = 0, i = 0; i < m; i++) {
			indexpos += indexedrows_rs[i].rank(startx);
			leaves[i] = cst.select_leaf(cst.csa.isa[indexpos] + 1);
			if (indexedrows_rs[i].rank(startx) != 0) {
				leavesmap[cst.lb(leaves[i])] = i;
			}
			indexpos += indexedrows_rs[i].rank(n) - indexedrows_rs[i].rank(startx) + 1;
		}

		for (size_type x = startx; x <= endx; x++) {
			/*std::cerr << "Computing f[x], with x = " << x << std::endl << std::flush;

			  std::cerr << "leaves for f[" << x << "] are: ";
			  for (auto l : leaves)
			  std::cerr << l << " ";
			  std::cerr << std::endl;*/

			size_type fimax = x;

			// Process each set of contiguous leaves
			for (size_type i = 0; i < m; i++) {
				node_t const l = leaves[i];
				if (!disable_efg_tricks && indexedrows_rs[i].rank(x) == 0)
					continue;
				if (cst.lb(l) == 0 || leavesmap.find(cst.lb(l) - 1) == leavesmap.end()) {
					// if leftmost leaf does not correspond to row i, skip
					if (concatenated_rows_rs.rank(cst.sn(cst.select_leaf(cst.lb(l) + 1))) != i)
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
								size_type gg = indexedrows_rs[ii].rank(x) + g;
								size_type fi;
								if (gg > indexedrows_rs[ii].rank(n)) {
									if (!disable_efg_tricks)
										fi = indexedrows_ss[ii].select(indexedrows_rs[ii].rank(n));
										// fi can be less than x here but it's still correct
									else
										fi = n;
									// fi can be less than x here but it's still correct
								} else {
									fi = indexedrows_ss[ii].select(gg);
								}
								// filter for first occurrence of ignore char
								if (ignorechars.length() > 0 && rowsignore_rs[ii].rank(x) != rowsignore_rs[ii].rank(n))
									fi = std::min(rowsignore_ss[ii].select(rowsignore_rs[ii].rank(x) + 1), fi);
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
			f[x] = std::max(f[x], fimax);

			for (size_type i = 0; i < m; i++) {
				if (MSA[i][x] != '-') {
					leavesmap.erase(cst.lb(leaves[i]));
					leaves[i] = cst.sl(leaves[i]);
					leavesmap[cst.lb(leaves[i])] = i;
				}
			}
		}
	}

	void compute_f(
			size_type const m,
			size_type const n,
			std::vector<size_type> &f, // store result here
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::string const &ignorechars,
			std::vector<sdsl::rank_support_v5<>> const &rowsignore_rs,
			std::vector<sdsl::select_support_mcl<>> const &rowsignore_ss,
			sdsl::rank_support_v5<> const &concatenated_rows_rs,
			std::vector<sdsl::rank_support_v5<>> const &indexedrows_rs,
			std::vector<sdsl::select_support_mcl<>> const &indexedrows_ss,
			const bool disable_efg_tricks = false
			) {
		/* Find the nodes corresponding to reading each whole row */
		std::vector<node_t> leaves(m, cst.root());
		std::unordered_map<size_type, size_type>  leavesmap;
		for (size_type next = 0, i = 0; i < m; i++) {
			leaves[i] = cst.select_leaf(cst.csa.isa[next] + 1);
			leavesmap[cst.lb(leaves[i])] = i;
			next += indexedrows_rs[i].rank(n) + 1;
		}

		/* binary coloring of the leaves */
		sdsl::bit_vector color(cst.size(cst.root()), false);
		sdsl::bit_vector fullrow;
		if (disable_efg_tricks)
			fullrow = sdsl::bit_vector(m, false); // mark if leaves[i] is still the initial value
		else
			fullrow = sdsl::bit_vector(m, true); // mark if leaves[i] is still the initial value

		for (size_type x = 0; x < n; x++) {
			/*std::cerr << "Computing f[x], with x = " << x << std::endl << std::flush;

			  std::cerr << "leaves for f[" << x << "] are: ";
			  for (auto l : leaves)
			  std::cerr << l << " ";
			  std::cerr << std::endl;*/

			size_type fimax = x;
			// Mark each leaf in leaves, filtering first all leaves corresponding to full rows
			for (size_type i = 0; i < m; i++) {
				if (fullrow[i])
					continue;
				for (size_type ll = cst.lb(leaves[i]); ll <= cst.rb(leaves[i]); ll++) { // cst is not a generalized suffix tree
					color[ll] = true;
				}
			}

			// Process each set of contiguous leaves
			for (size_type i = 0; i < m; i++) {
				node_t const l = leaves[i];
				if (fullrow[i])
					continue;
				if (cst.lb(l) == 0 || color[cst.lb(l) - 1] == false) {
					// if leftmost leaf does not correspond to row i, skip
					if (concatenated_rows_rs.rank(cst.sn(cst.select_leaf(cst.lb(l) + 1))) != i)
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
								size_type gg = indexedrows_rs[ii].rank(x) + g;
								size_type fi;
								if (gg > indexedrows_rs[ii].rank(n)) {
									if (!disable_efg_tricks)
										fi = indexedrows_ss[ii].select(indexedrows_rs[ii].rank(n));
										// fi can be less than x here but it's still correct
									else
										fi = n;
								} else {
									fi = indexedrows_ss[ii].select(gg);
								}
								// filter for first occurrence of ignore char
								if (ignorechars.length() > 0 && rowsignore_rs[ii].rank(x) != rowsignore_rs[ii].rank(n))
									fi = std::min(rowsignore_ss[ii].select(rowsignore_rs[ii].rank(x) + 1), fi);
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
			f[x] = std::max(f[x], fimax);

			for (size_type i = 0; i < m; i++) {
				for (size_type ll = cst.lb(leaves[i]); ll <= cst.rb(leaves[i]); ll++) { // cst is not a generalized suffix tree
					color[ll] = false;
				}
				if (MSA[i][x] != '-') {
					leavesmap.erase(cst.lb(leaves[i]));
					leaves[i] = cst.sl(leaves[i]);
					leavesmap[cst.lb(leaves[i])] = i;
					fullrow[i] = false;
				}
			}
		}
	}

	size_type count_interval_union(
		const size_type m,
		const std::vector<size_type> &l,
		const std::vector<size_type> &r,
		const std::vector<bool> &to_consider,
		const std::vector<bool> &to_ignore
		)
	{
		std::vector<std::pair<size_type,size_type>> intervals;
		for (size_type i = 0; i < m; i++) {
			if (to_consider[i] && !to_ignore[i])
				intervals.push_back({ l[i], r[i] });
		}

		// sort by starting interval pos
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<0>(a) < std::get<0>(b));
			}
		};
		std::sort(intervals.begin(), intervals.end(), Local::pair_comparator);

		size_type count = 0;
		size_type active_interval_end_pos = 0; // first position where active interval is no longer active
		for (const auto [start, end] : intervals) {
			if (start >= active_interval_end_pos) {
				count += end - start + 1;
				active_interval_end_pos = end + 1;
			} else if (end < active_interval_end_pos) { // && start < active_interval_end_pos
				// fully contained in active interval, do nothing
			} else { // start < active_interval_end_pos && active_interval_end_pos <= end
				// intersects with active interval
				count += end - active_interval_end_pos + 1;
				active_interval_end_pos = end + 1;
			}
		}

		return count;
	}

	void compute_f_heuristic(
			size_type const m,
			size_type const n,
			size_type const x,
			std::vector<size_type> &f, // store result here
			std::vector<std::string> const &MSA,
			sa_type const &sa,
			std::string const &ignorechars,
			std::vector<sdsl::rank_support_v5<>> const &rowsignore_rs,
			sdsl::rank_support_v5<> const &concatenated_rows_rs,
			std::vector<sdsl::rank_support_v5<>> const &indexedrows_rs,
			std::vector<sdsl::select_support_mcl<>> const &indexedrows_ss,
			const bool disable_efg_tricks = false
			) {
		//std::cerr << "Call to compute_f_heuristic for x = " << x << ", current f[x] is " << f[x] << std::endl;
		// find ranges corresponding to strings starting at the (x+1)-th col and ending at the f[x]-th col
		std::vector<size_type> l(m, 0);
		std::vector<size_type> r(m, sa.size() - 1);
		std::vector<bool> initialized(m, false);
		std::vector<bool> to_ignore(m, false);
		size_type active_rows = 0;
		for (size_type i = 0; i < m; i++) {
			if (indexedrows_rs[i].rank(x) != 0) {
				initialized[i] = true;
				active_rows += 1;
				if (indexedrows_rs[i].rank(x) == indexedrows_rs[i].rank(n)) {
				} else {
					f[x] = std::max(f[x], indexedrows_ss[i](indexedrows_rs[i].rank(x)+1));
				}
			} else {
				initialized[i] = false;
			}
		}
		for (size_type i = 0; i < m; i++) {
			if (initialized[i]) {
				//std::cerr << "row " << i << "\nlooking for string MSA[" << i << "," << x << ".." << f[x] << "\n";
				/*for (size_type j = x; j <= f[x]; j++) {
					if (MSA[i][j] != '-')
						std::cerr << MSA[i][j];
				}
				std::cerr << std::endl;*/

				std::string s = remove_gaps(MSA[i].substr(x,f[x] - x + 1));
				int res = sdsl::forward_search(sa, l[i], r[i], s.begin(), s.end(), l[i], r[i]);
				assert(res != 0);
			}
		}

		int iterations = 0;
		while (f[x] < n - 1 && count_interval_union(m, l, r, initialized, to_ignore) > active_rows) {
			f[x] += f[x] - x + 1;
			iterations += 1;
			if (iterations >= 5 || f[x] >= n - 1 || f[x] - x >= 50000) {// TODO disable efg tricks here?
				f[x] = n - 1;
				break;
			}
			for (size_type i = 0; i < m; i++) {
				if (!to_ignore[i] && MSA[i][f[x]] != '-') {
					if (!initialized[i]) {
						active_rows += 1;
						initialized[i] = true;
					}
					if (ignorechars.size() > 0 && rowsignore_rs[i](f[x]) - rowsignore_rs[i](f[x]-1) > 0) {
						to_ignore[i] = true;
						if (initialized[i]) {
							active_rows -= 1; // we lose optimality!?
						}
					} else {
						std::string s = remove_gaps(MSA[i].substr(x,f[x] - x + 1));
						int res = sdsl::forward_search(sa, l[i], r[i], s.begin(), s.end(), l[i], r[i]);
						assert(res != 0);
					}
				}
			}
		}
	}

	void compute_f_heuristic_interleaved(
			size_type const m,
			size_type const n,
			size_type const startx,
			size_type const jump,
			std::vector<size_type> &f, // store result here
			std::vector<std::string> const &MSA,
			sa_type const &sa,
			std::string const &ignorechars,
			std::vector<sdsl::rank_support_v5<>> const &rowsignore_rs,
			sdsl::rank_support_v5<> const &concatenated_rows_rs,
			std::vector<sdsl::rank_support_v5<>> const &indexedrows_rs,
			std::vector<sdsl::select_support_mcl<>> const &indexedrows_ss,
			const bool disable_efg_tricks = false
			) {
		//std::cerr << "Call to compute_f_heuristic for x = " << x << ", current f[x] is " << f[x] << std::endl;
		// find ranges corresponding to strings starting at the (x+1)-th col and ending at the f[x]-th col
		for (size_type x = startx; x < n; x += jump) {
			compute_f_heuristic(m, n, x, f, MSA, sa, ignorechars, rowsignore_rs, concatenated_rows_rs, indexedrows_rs, indexedrows_ss, disable_efg_tricks);
		}
	}

	void segment_elastic_minmaxlength(
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::string &ignorechars,
			std::vector<size_type> &out_indices,
			const bool disable_efg_tricks,
			std::vector<size_type> &f,
			const bool segment = true
			) {
		size_type const n = MSA[0].size();
		size_type const m = MSA.size();
		size_type nongap = 0;
		size_type toignore = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					nongap += 1;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					toignore += 1;
			}
		}

		std::cerr << "MSA contains " << (n*m) - nongap  << " gaps.\n" << std::flush;
		std::cerr << "MSA contains " << toignore  << " characters to ignore for the semi-repeat-free property.\n" << std::flush;

		// Preprocess the MSA rows for rank, select queries
		std::vector<sdsl::bit_vector> indexedrows(m, sdsl::bit_vector(n, 1));
		std::vector<sdsl::bit_vector> rowsignore(m, sdsl::bit_vector(n, 0)); // ignore chars
		for (size_type i = 0; i < m; i++) {
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] == '-')
					indexedrows[i][j] = 0;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					rowsignore[i][j] = 1;
			}
		}

		// rank support on the rows
		std::vector<sdsl::rank_support_v5<>> indexedrows_rs;
		indexedrows_rs.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::rank_support_v5<> rs(&indexedrows[i]);
			indexedrows_rs[i] = rs;
		}

		// select support on the rows
		std::vector<sdsl::select_support_mcl<>> indexedrows_ss;
		indexedrows_ss.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::select_support_mcl<> ss(&indexedrows[i]);
			indexedrows_ss[i] = ss;
		}

		// rank, select on ignore chars
		std::vector<sdsl::rank_support_v5<>> rowsignore_rs;
		if (ignorechars.length() > 0) {
			rowsignore_rs.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::rank_support_v5<> rs(&rowsignore[i]);
				rowsignore_rs[i] = rs;
			}
		}
		std::vector<sdsl::select_support_mcl<>> rowsignore_ss;
		if (ignorechars.length() > 0) {
			rowsignore_ss.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::select_support_mcl<> ss(&rowsignore[i]);
				rowsignore_ss[i] = ss;
			}
		}

		// bitvector + rank and select support for finding corresponding row
		sdsl::bit_vector concatenated_rows(nongap + m, 0);
		size_type k = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					k += 1;
			}
			concatenated_rows[k++] = 1;
		}
		sdsl::rank_support_v5<> concatenated_rows_rs(&concatenated_rows);

		/* This is a preprocessing part of the main algorithm */

		/* f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is a semi-repeat-free segment */
		//f = std::vector<size_type>(n, 0);
		compute_f(m, n, std::ref(f), std::ref(MSA), std::ref(cst), std::ref(ignorechars), std::ref(rowsignore_rs), std::ref(rowsignore_ss), std::ref(concatenated_rows_rs), std::ref(indexedrows_rs), std::ref(indexedrows_ss), disable_efg_tricks);

		/*std::cerr << "f: ";
		  for (auto fx : f)
		  std::cerr << fx << " ";
		  std::cerr << std::endl;*/
		if (!segment)
			return;

		if (disable_efg_tricks) {
			if (f[0] == n) {
				std::cerr << "No valid segmentation found!\n";
				exit(1);
			}
		}
		// else there is always a valid segmentation

		// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
		std::vector<std::pair<size_type,size_type>> minimal_right_extensions;
		minimal_right_extensions.resize(n);
		for (size_type x = 0; x < n; x++) {
			std::pair<size_type,size_type> p(x, f[x]+1);
			minimal_right_extensions[x] = p;
		}
		// TODO: choose sorting algorithm
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<1>(a) < std::get<1>(b));
			}
		};
		std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), Local::pair_comparator);


		// END OF PREPROCESSING

		std::cerr << "Computing optimal segmentation..." << std::flush;
		// TODO: swap size_type with int32 or optimal multiple of 2
		std::vector<size_type> count_solutions(n, 0);
		std::vector<size_type> backtrack_count(n, 0);
		std::vector<std::list<std::pair<size_type,size_type>>> transition_list(n + 2);
		// NOTE: transition_list[n + 1] can be used to manage the terminator character
		// UPDATE: this note probably is not true, this can be fixed in here
		std::vector<size_type> minmaxlength(n + 1, 0);
		std::vector<size_type> backtrack(n + 1, 0);
		size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
		for (size_type j = 1; j <= n; j++) {
			while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
				size_type xy  = std::get<0>(minimal_right_extensions[y]);
				size_type rec_score = minmaxlength[xy];
				//std::cerr << "rec_score = " << rec_score << std::endl;
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
			if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S? should we pick the smallest backtrack?
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
		/*std::cerr << "minmaxlength: ";
		  for (auto v : minmaxlength)
		  std::cerr << v << " ";
		  std::cerr << std::endl;
		  std::cerr << "backtrack: ";
		  for (auto v : backtrack)
		  std::cerr << v << " ";
		  std::cerr << std::endl;*/
		std::cerr << "done (optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;

		// NB: added block info
		std::list<size_type> boundariestemp;
		size_type j = n;
		boundariestemp.push_front(j);
		while (backtrack[j]!=0) {
			//std::cerr << j << " ";
			boundariestemp.push_front(backtrack[j]-1);
			j = backtrack[j];
		}

		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);

		std::swap(boundaries, out_indices);
	}

	void segment_elastic_minmaxlength(
			std::vector<size_type> &out_indices,
			const bool disable_efg_tricks,
			std::vector<size_type> &f,
			size_type n
			) {
		if (disable_efg_tricks) {
			if (f[0] == n) {
				std::cerr << "No valid segmentation found!\n";
				exit(1);
			}
		}
		// else there is always a valid segmentation

		// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
		std::vector<std::pair<size_type,size_type>> minimal_right_extensions;
		minimal_right_extensions.resize(n);
		for (size_type x = 0; x < n; x++) {
			std::pair<size_type,size_type> p(x, f[x]+1);
			minimal_right_extensions[x] = p;
		}
		// TODO: choose sorting algorithm
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<1>(a) < std::get<1>(b));
			}
		};
		std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), Local::pair_comparator);


		// END OF PREPROCESSING

		std::cerr << "Computing optimal segmentation..." << std::flush;
		// TODO: swap size_type with int32 or optimal multiple of 2
		std::vector<size_type> count_solutions(n, 0);
		std::vector<size_type> backtrack_count(n, 0);
		std::vector<std::list<std::pair<size_type,size_type>>> transition_list(n + 2);
		// NOTE: transition_list[n + 1] can be used to manage the terminator character
		// UPDATE: this note probably is not true, this can be fixed in here
		std::vector<size_type> minmaxlength(n + 1, 0);
		std::vector<size_type> backtrack(n + 1, 0);
		size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
		for (size_type j = 1; j <= n; j++) {
			while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
				size_type xy  = std::get<0>(minimal_right_extensions[y]);
				size_type rec_score = minmaxlength[xy];
				//std::cerr << "rec_score = " << rec_score << std::endl;
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
			if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S? should we pick the smallest backtrack?
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
		/*std::cerr << "minmaxlength: ";
		  for (auto v : minmaxlength)
		  std::cerr << v << " ";
		  std::cerr << std::endl;
		  std::cerr << "backtrack: ";
		  for (auto v : backtrack)
		  std::cerr << v << " ";
		  std::cerr << std::endl;*/
		std::cerr << "done (optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;

		// NB: added block info
		std::list<size_type> boundariestemp;
		size_type j = n;
		boundariestemp.push_front(j);
		while (backtrack[j]!=0) {
			//std::cerr << j << " ";
			boundariestemp.push_front(backtrack[j]-1);
			j = backtrack[j];
		}

		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);

		std::swap(boundaries, out_indices);
	}

	void segment_elastic_minmaxlength_worker(
			char const *input_path,
			std::string &ignorechars,
			std::vector<size_type> &out_indices,
			const bool disable_efg_tricks,
			std::vector<size_type> &f,
			const long rows,
			const size_type m,
			const size_type n
			) {
		cst_type cst;
		int miniMSArow = -1;
		std::vector<std::string> miniMSA;
		miniMSA = load_rows(input_path, rows, miniMSArow);
		while (miniMSA.size() > 0) {
			compute_cst_im(miniMSA, cst);
			segment_elastic_minmaxlength(miniMSA, cst, ignorechars, out_indices, disable_efg_tricks, f, false);

			miniMSA = load_rows(input_path, rows, miniMSArow);
		}
	}

	void segment_elastic_minmaxlength_multithread(
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::string &ignorechars,
			std::vector<size_type> &out_indices,
			int max_threads,
			const bool disable_efg_tricks,
			std::vector<size_type> &f
			) {
		size_type const n = MSA[0].size();
		size_type const m = MSA.size();
		size_type nongap = 0;
		size_type toignore = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					nongap += 1;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					toignore += 1;
			}
		}

		std::cerr << "MSA contains " << (n*m) - nongap  << " gaps.\n" << std::flush;
		std::cerr << "MSA contains " << toignore  << " characters to ignore for the semi-repeat-free property.\n" << std::flush;

		// Preprocess the MSA rows for rank, select queries
		std::vector<sdsl::bit_vector> indexedrows(m, sdsl::bit_vector(n, 1));
		std::vector<sdsl::bit_vector> rowsignore(m, sdsl::bit_vector(n, 0)); // ignore chars
		for (size_type i = 0; i < m; i++) {
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] == '-')
					indexedrows[i][j] = 0;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					rowsignore[i][j] = 1;
			}
		}

		// rank support on the rows
		std::vector<sdsl::rank_support_v5<>> indexedrows_rs;
		indexedrows_rs.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::rank_support_v5<> rs(&indexedrows[i]);
			indexedrows_rs[i] = rs;
		}

		// select support on the rows
		std::vector<sdsl::select_support_mcl<>> indexedrows_ss;
		indexedrows_ss.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::select_support_mcl<> ss(&indexedrows[i]);
			indexedrows_ss[i] = ss;
		}

		// rank, select on ignore chars
		std::vector<sdsl::rank_support_v5<>> rowsignore_rs;
		if (ignorechars.length() > 0) {
			rowsignore_rs.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::rank_support_v5<> rs(&rowsignore[i]);
				rowsignore_rs[i] = rs;
			}
		}
		std::vector<sdsl::select_support_mcl<>> rowsignore_ss;
		if (ignorechars.length() > 0) {
			rowsignore_ss.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::select_support_mcl<> ss(&rowsignore[i]);
				rowsignore_ss[i] = ss;
			}
		}

		// bitvector + rank and select support for finding corresponding row
		sdsl::bit_vector concatenated_rows(nongap + m, 0);
		size_type k = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					k += 1;
			}
			concatenated_rows[k++] = 1;
		}
		sdsl::rank_support_v5<> concatenated_rows_rs(&concatenated_rows);

		/* This is a preprocessing part of the main algorithm */
		/* Find the nodes corresponding to reading each whole row */
		std::vector<node_t> leaves(m, cst.root());
		std::unordered_map<size_type, size_type>  leavesmap;
		for (size_type next = 0, i = 0; i < m; i++) {
			leaves[i] = cst.select_leaf(cst.csa.isa[next] + 1);
			leavesmap[cst.lb(leaves[i])] = i;
			next += indexedrows_rs[i].rank(n) + 1;
		}

		/* binary coloring of the leaves */
		sdsl::bit_vector color(cst.size(cst.root()), false);

		/* f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is a semi-repeat-free segment */
		//f = std::vector<size_type>(n, 0);
		for (size_type x = 0; x < n;) {
			std::vector<std::thread> t;
			for (int i = 0; i < max_threads && x < n; i++) {
				//std::cerr << "Calling compute_f_range for range" << x << " to " << std::min(x + n/max_threads, n-1) << std::endl << std::flush;
				t.push_back(std::thread(compute_f_range, m, n, x, std::min(x + n/max_threads, n-1), std::ref(f), std::ref(MSA), std::ref(cst), std::ref(ignorechars), std::ref(rowsignore_rs), std::ref(rowsignore_ss), std::ref(concatenated_rows_rs), std::ref(indexedrows_rs), std::ref(indexedrows_ss), disable_efg_tricks));
				x = std::min(x + n/max_threads, n-1) + 1;
			}
			for (auto &tt : t)
			{
				tt.join();
			}
		}

		if (disable_efg_tricks) {
			if (f[0] == n) {
				std::cerr << "No valid segmentation found!\n";
				exit(1);
			}
		}
		// else there is always a valid segmentation

		// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
		std::vector<std::pair<size_type,size_type>> minimal_right_extensions;
		minimal_right_extensions.resize(n);
		for (size_type x = 0; x < n; x++) {
			std::pair<size_type,size_type> p(x, f[x]+1);
			minimal_right_extensions[x] = p;
		}
		// TODO: choose sorting algorithm
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<1>(a) < std::get<1>(b));
			}
		};
		std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), Local::pair_comparator);


		// END OF PREPROCESSING

		std::cerr << "Computing optimal segmentation..." << std::flush;
		std::vector<size_type> count_solutions(n, 0);
		std::vector<size_type> backtrack_count(n, 0);
		std::vector<std::list<std::pair<size_type,size_type>>> transition_list(n + 2);
		// NOTE: transition_list[n + 1] can be used to manage the terminator character
		// UPDATE: this note probably is not true, this can be fixed in here
		std::vector<size_type> minmaxlength(n + 1, 0);
		std::vector<size_type> backtrack(n + 1, 0);
		size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
		for (size_type j = 1; j <= n; j++) {
			while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
				size_type xy  = std::get<0>(minimal_right_extensions[y]);
				size_type rec_score = minmaxlength[xy];
				//std::cerr << "rec_score = " << rec_score << std::endl;
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
			if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S? should we pick the smallest backtrack?
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
		/*std::cerr << "minmaxlength: ";
		  for (auto v : minmaxlength)
		  std::cerr << v << " ";
		  std::cerr << std::endl;
		  std::cerr << "backtrack: ";
		  for (auto v : backtrack)
		  std::cerr << v << " ";
		  std::cerr << std::endl;*/
		std::cerr << "done (optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;

		// NB: added block info
		std::list<size_type> boundariestemp;
		size_type j = n;
		boundariestemp.push_front(j);
		while (backtrack[j]!=0) {
			//std::cerr << j << " ";
			boundariestemp.push_front(backtrack[j]-1);
			j = backtrack[j];
		}

		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);

		std::swap(boundaries, out_indices);
	}

	void segment_elastic_minmaxlength_heuristic(
			std::vector<std::string> const &MSA,
			sa_type const &sa,
			std::string &ignorechars,
			std::vector<size_type> &out_indices,
			const bool disable_efg_tricks,
			std::vector<size_type> &f
			) {
		size_type const n = MSA[0].size();
		size_type const m = MSA.size();
		size_type nongap = 0;
		size_type toignore = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					nongap += 1;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					toignore += 1;
			}
		}

		std::cerr << "MSA contains " << (n*m) - nongap  << " gaps.\n" << std::flush;
		std::cerr << "MSA contains " << toignore  << " characters to ignore for the semi-repeat-free property.\n" << std::flush;

		// Preprocess the MSA rows for rank, select queries
		std::cerr << "preprocessing MSA...";
		std::vector<sdsl::bit_vector> indexedrows(m, sdsl::bit_vector(n, 1));
		std::vector<sdsl::bit_vector> rowsignore(m, sdsl::bit_vector(n, 0)); // ignore chars
		for (size_type i = 0; i < m; i++) {
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] == '-')
					indexedrows[i][j] = 0;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					rowsignore[i][j] = 1;
			}
		}

		// rank support on the rows
		std::vector<sdsl::rank_support_v5<>> indexedrows_rs;
		indexedrows_rs.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::rank_support_v5<> rs(&indexedrows[i]);
			indexedrows_rs[i] = rs;
		}

		// select support on the rows
		std::vector<sdsl::select_support_mcl<>> indexedrows_ss;
		indexedrows_ss.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::select_support_mcl<> ss(&indexedrows[i]);
			indexedrows_ss[i] = ss;
		}

		// rank, select on ignore chars
		std::vector<sdsl::rank_support_v5<>> rowsignore_rs;
		if (ignorechars.length() > 0) {
			rowsignore_rs.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::rank_support_v5<> rs(&rowsignore[i]);
				rowsignore_rs[i] = rs;
			}
		}

		// bitvector + rank and select support for finding corresponding row
		sdsl::bit_vector concatenated_rows(nongap + m, 0);
		size_type k = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					k += 1;
			}
			concatenated_rows[k++] = 1;
		}
		sdsl::rank_support_v5<> concatenated_rows_rs(&concatenated_rows);
		std::cerr << "done." << std::endl;

		/* This is a preprocessing part of the main algorithm */

		/* f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is a semi-repeat-free segment */
		std::cerr << "Computing f values...";
		for (size_type x = 0; x < n; x++) {
			compute_f_heuristic(m, n, x, std::ref(f), std::ref(MSA), std::ref(sa), std::ref(ignorechars), std::ref(rowsignore_rs), std::ref(concatenated_rows_rs), std::ref(indexedrows_rs), std::ref(indexedrows_ss), disable_efg_tricks);
		}
		std::cerr << "done." << std::endl;

		/*std::cerr << "f: ";
		  for (auto fx : f)
		  std::cerr << fx << " ";
		  std::cerr << std::endl;*/

		if (disable_efg_tricks) {
			if (f[0] == n) {
				std::cerr << "No valid segmentation found!\n";
				exit(1);
			}
		}
		// else there is always a valid segmentation

		std::cerr << "Computing optimal segmentation...";
		// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
		std::vector<std::pair<size_type,size_type>> minimal_right_extensions;
		minimal_right_extensions.resize(n);
		for (size_type x = 0; x < n; x++) {
			std::pair<size_type,size_type> p(x, f[x]+1);
			minimal_right_extensions[x] = p;
		}
		// TODO: choose sorting algorithm
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<1>(a) < std::get<1>(b));
			}
		};
		std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), Local::pair_comparator);


		// END OF PREPROCESSING

		// TODO: swap size_type with int32 or optimal multiple of 2
		std::vector<size_type> count_solutions(n, 0);
		std::vector<size_type> backtrack_count(n, 0);
		std::vector<std::list<std::pair<size_type,size_type>>> transition_list(n + 2);
		// NOTE: transition_list[n + 1] can be used to manage the terminator character
		// UPDATE: this note probably is not true, this can be fixed in here
		std::vector<size_type> minmaxlength(n + 1, 0);
		std::vector<size_type> backtrack(n + 1, 0);
		size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
		for (size_type j = 1; j <= n; j++) {
			while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
				size_type xy  = std::get<0>(minimal_right_extensions[y]);
				size_type rec_score = minmaxlength[xy];
				//std::cerr << "rec_score = " << rec_score << std::endl;
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
			if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S? should we pick the smallest backtrack?
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
		/*std::cerr << "minmaxlength: ";
		  for (auto v : minmaxlength)
		  std::cerr << v << " ";
		  std::cerr << std::endl;
		  std::cerr << "backtrack: ";
		  for (auto v : backtrack)
		  std::cerr << v << " ";
		  std::cerr << std::endl;*/
		std::cerr << "done (optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;

		// NB: added block info
		std::list<size_type> boundariestemp;
		size_type j = n;
		boundariestemp.push_front(j);
		while (backtrack[j]!=0) {
			//std::cerr << j << " ";
			boundariestemp.push_front(backtrack[j]-1);
			j = backtrack[j];
		}

		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);

		std::swap(boundaries, out_indices);
	}

	void segment_elastic_minmaxlength_heuristic_multithread(
			std::vector<std::string> const &MSA,
			sa_type const &sa,
			std::string &ignorechars,
			std::vector<size_type> &out_indices,
			int max_threads,
			const bool disable_efg_tricks,
			std::vector<size_type> &f
			) {
		size_type const n = MSA[0].size();
		size_type const m = MSA.size();
		size_type nongap = 0;
		size_type toignore = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					nongap += 1;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					toignore += 1;
			}
		}

		std::cerr << "MSA contains " << (n*m) - nongap  << " gaps.\n" << std::flush;
		std::cerr << "MSA contains " << toignore  << " characters to ignore for the semi-repeat-free property.\n" << std::flush;

		// Preprocess the MSA rows for rank, select queries
		std::cerr << "preprocessing MSA...";
		std::vector<sdsl::bit_vector> indexedrows(m, sdsl::bit_vector(n, 1));
		std::vector<sdsl::bit_vector> rowsignore(m, sdsl::bit_vector(n, 0)); // ignore chars
		for (size_type i = 0; i < m; i++) {
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] == '-')
					indexedrows[i][j] = 0;
				if (ignorechars.find(MSA[i][j]) != std::string::npos)
					rowsignore[i][j] = 1;
			}
		}

		// rank support on the rows
		std::vector<sdsl::rank_support_v5<>> indexedrows_rs;
		indexedrows_rs.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::rank_support_v5<> rs(&indexedrows[i]);
			indexedrows_rs[i] = rs;
		}

		// select support on the rows
		std::vector<sdsl::select_support_mcl<>> indexedrows_ss;
		indexedrows_ss.resize(m);
		for (size_type i = 0; i < m; i++) {
			sdsl::select_support_mcl<> ss(&indexedrows[i]);
			indexedrows_ss[i] = ss;
		}

		// rank, select on ignore chars
		std::vector<sdsl::rank_support_v5<>> rowsignore_rs;
		if (ignorechars.length() > 0) {
			rowsignore_rs.resize(m);
			for (size_type i = 0; i < m; i++) {
				sdsl::rank_support_v5<> rs(&rowsignore[i]);
				rowsignore_rs[i] = rs;
			}
		}

		// bitvector + rank and select support for finding corresponding row
		sdsl::bit_vector concatenated_rows(nongap + m, 0);
		size_type k = 0;
		for (size_type i = 0; i < m; i++) { // TODO: MSA struct
			for (size_type j = 0; j < n; j++) {
				if (MSA[i][j] != '-')
					k += 1;
			}
			concatenated_rows[k++] = 1;
		}
		sdsl::rank_support_v5<> concatenated_rows_rs(&concatenated_rows);
		std::cerr << "done." << std::endl;

		/* This is a preprocessing part of the main algorithm */

		/* f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is a semi-repeat-free segment */
		std::cerr << "Computing f values...";
		for (size_type x = 0; x < n;) {
			std::vector<std::thread> t;
			for (int i = 0; i < max_threads && x < n; i++) {
				//std::cerr << "Calling compute_f_range for range" << x << " to " << std::min(x + n/max_threads, n-1) << std::endl << std::flush;
				t.push_back(std::thread(compute_f_heuristic_interleaved, m, n, i, max_threads, std::ref(f), std::ref(MSA), std::ref(sa), std::ref(ignorechars), std::ref(rowsignore_rs), std::ref(concatenated_rows_rs), std::ref(indexedrows_rs), std::ref(indexedrows_ss), disable_efg_tricks));
			}
			for (auto &tt : t)
			{
				tt.join();
			}
		}
		std::cerr << "done." << std::endl;

		/*std::cerr << "f: ";
		  for (auto fx : f)
		  std::cerr << fx << " ";
		  std::cerr << std::endl;*/

		if (disable_efg_tricks) {
			if (f[0] == n) {
				std::cerr << "No valid segmentation found!\n";
				exit(1);
			}
		}
		// else there is always a valid segmentation

		std::cerr << "Computing optimal segmentation..." << std::flush;
		// Sort the resulting pairs (x,f(x)), make f(x) 1-indexed
		std::vector<std::pair<size_type,size_type>> minimal_right_extensions;
		minimal_right_extensions.resize(n);
		for (size_type x = 0; x < n; x++) {
			std::pair<size_type,size_type> p(x, f[x]+1);
			minimal_right_extensions[x] = p;
		}
		// TODO: choose sorting algorithm
		struct Local {
			static bool pair_comparator(std::pair<size_type,size_type> a, std::pair<size_type,size_type> b) {
				return (std::get<1>(a) < std::get<1>(b));
			}
		};
		std::sort(minimal_right_extensions.begin(), minimal_right_extensions.end(), Local::pair_comparator);


		// END OF PREPROCESSING

		// TODO: swap size_type with int32 or optimal multiple of 2
		std::vector<size_type> count_solutions(n, 0);
		std::vector<size_type> backtrack_count(n, 0);
		std::vector<std::list<std::pair<size_type,size_type>>> transition_list(n + 2);
		// NOTE: transition_list[n + 1] can be used to manage the terminator character
		// UPDATE: this note probably is not true, this can be fixed in here
		std::vector<size_type> minmaxlength(n + 1, 0);
		std::vector<size_type> backtrack(n + 1, 0);
		size_type y = 0, I = 0, S = n + 1, backtrack_S = -1;
		for (size_type j = 1; j <= n; j++) {
			while (y < n && j == std::get<1>(minimal_right_extensions[y])) {
				size_type xy  = std::get<0>(minimal_right_extensions[y]);
				size_type rec_score = minmaxlength[xy];
				//std::cerr << "rec_score = " << rec_score << std::endl;
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
			if (count_solutions[I] > 0 && I < S) { //TODO: what if I == S? should we pick the smallest backtrack?
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
		/*std::cerr << "minmaxlength: ";
		  for (auto v : minmaxlength)
		  std::cerr << v << " ";
		  std::cerr << std::endl;
		  std::cerr << "backtrack: ";
		  for (auto v : backtrack)
		  std::cerr << v << " ";
		  std::cerr << std::endl;*/
		std::cerr << "done (optimal segment length = " << minmaxlength[n] << ")." << std::endl << std::flush;

		// NB: added block info
		std::list<size_type> boundariestemp;
		size_type j = n;
		boundariestemp.push_front(j);
		while (backtrack[j]!=0) {
			//std::cerr << j << " ";
			boundariestemp.push_front(backtrack[j]-1);
			j = backtrack[j];
		}

		std::vector<size_type> boundaries;
		for (const auto& j : boundariestemp)
			boundaries.push_back(j);

		std::swap(boundaries, out_indices);
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
		std::string buffer;
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

				buffer.clear();
				buffer += src_label;
				buffer += dst_label;
				buffer += fbg::g_separator_character;
				std::reverse(buffer.begin(), buffer.end());

				if (VERBOSE_LOGGING)
					std::cerr << "Outputting “" << buffer << "”\n";

				temp_os << buffer;

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
		// These are only used when debugging.
		std::map <size_type, std::string> labels_b;
		std::map <size_type, std::string> labels_e;
		for (std::size_t i(0), count(node_labels.size()); i < count; ++i)
		{
			auto const &label(node_labels[i]);
			if (VERBOSE_LOGGING)
				std::cerr << "Handling node label “" << label << "”\n";

			// The search direction is forward as the index was generated from reversed labels.
			fbg::csa_type::size_type lhs{}, rhs{};
			auto const match_count(
					sdsl::backward_search(
						csa,
						0,
						csa_size - 1,
						label.crbegin(),
						label.crend(),
						lhs,
						rhs
						)
					);

			assert(match_count);
			assert(lhs < b_positions.size());
			assert(rhs < e_positions.size());

			// Sanity check.
			{
				bool should_stop(false);

				if (b_positions[lhs])
				{
					should_stop = true;
					std::cerr << "b_positions[lhs] already set for " << lhs << ", edge = " << label << '\n';
					if (SHOULD_STORE_ASSIGNED_NODE_LABELS)
						std::cerr << "Previously assigned: " << labels_b[lhs] << '\n';
				}

				if (e_positions[rhs])
				{
					should_stop = true;
					std::cerr << "e_positions[rhs] already set for " << rhs << ", edge = " << label << '\n';
					if (SHOULD_STORE_ASSIGNED_NODE_LABELS)
						std::cerr << "Previously assigned: " << labels_e[rhs] << '\n';
				}

				// SIKE
				//if (should_stop)
				//	std::abort();
			}

			if (SHOULD_STORE_ASSIGNED_NODE_LABELS)
			{
				labels_b[lhs] = label;
				labels_e[rhs] = label;
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

	void make_gfa(
			int const m,
			int const n,
			std::vector<std::string> const &identifiers,
			std::vector <std::string> &node_labels,
			adjacency_list &edges,
			std::vector <size_type> &node_blocks,
			std::vector <size_type> &block_indices,
			bool output_paths,
			std::vector<std::vector<size_type>> &paths,
			char const *dst_path_c,
			char const *mem_chart_path_c)
	{
		// Create a temporary file.
		std::string const dst_path(dst_path_c);
		std::fstream dst_os;
		dst_os.exceptions(std::fstream::failbit);
		dst_os.open(dst_path, std::ios_base::out);

		assert(edges.size() == node_labels.size());
		// MSA info
		dst_os << "M\t" << m << "\t" << n << std::endl;

		// Segmentation info
		dst_os << "X\t1";
		for (size_type i = 0; (int)i < (int)block_indices.size() - 1; i++)
			dst_os << "\t" << block_indices[i] + 2;
		dst_os << std::endl;

		// Partition of nodes into blocks
		dst_os << "B\t";
		size_t pastblock = -1, blockheight = 0;
		for (size_t i = 0; i < node_labels.size(); i++)
		{
			if (node_blocks[i] != pastblock && pastblock != (size_t)-1) {
				dst_os << blockheight << "\t";
				blockheight = 1;
			} else {
				blockheight += 1;
			}
			pastblock = node_blocks[i];
		}
		dst_os << blockheight << std::endl;

		// node labels
		for (size_t i = 0; i < node_labels.size(); ++i)
		{
			auto const &src_label(node_labels[i]);
			dst_os << "S\t" << i << "\t" << src_label << std::endl;
		}

		// edges
		for (size_t i = 0; i < node_labels.size(); i++)
		{
			// sort the destination nodes
			auto const &dst_nodes(edges[i]);
			std::vector<size_type> sorted_dst_nodes(dst_nodes.begin(), dst_nodes.end());
			std::sort(sorted_dst_nodes.begin(), sorted_dst_nodes.end());
			for (auto const dst_node : sorted_dst_nodes)
				dst_os << "L\t" << i << "\t+\t" << dst_node << "\t+\t0M" << std::endl;
		}

		if (!output_paths) return;
		// paths
		for (size_type i = 0; i < paths.size(); i++)
		{
			assert(m == paths.size());
			assert(identifiers.size() == paths.size());
			dst_os << "P\t" << identifiers[i] << "\t";
			for (size_type j = 0; j < paths[i].size() - 1; j++) {
				dst_os << paths[i][j] << "+,";
			}
			dst_os << paths[i][paths[i].size()-1] << "+";
			dst_os << "\t*" << std::endl;
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
				os << 'n' << i << " [label = \"" << i+1 << ": ";
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

	bool contains_ignore_chars(const std::string &s, const std::string &ignorechars)
	{
		for (char c : ignorechars) {
			if (s.find(c) != std::string::npos)
				return true;
		}
		return false;
	}


	// check semi-repeat-free property of a node
	bool efg_validate_node(
			const size_type node,
			const std::vector<std::string> &node_labels,
			const std::vector<std::pair<size_type,size_type>> &ordered_edges,
			const std::vector<size_type> &node_blocks,
			const std::vector<bool> &is_source,
			const std::vector<bool> &is_sink,
			const std::string &ignore_chars,
			const sdsl::csa_wt<> &index,
			const sdsl::rank_support_v5<> &dels_rs,
			const sdsl::select_support_mcl<> &dels_ss)
	{
		if (is_source[node] or is_sink[node])
			return true;

		if (ignore_chars.length() > 0 && contains_ignore_chars(node_labels[node], ignore_chars))
			return true;


		//auto occs = sdsl::locate(index, node_labels[node]);
		sdsl::csa_wt<>::size_type l, r;
		sdsl::backward_search(index, 0, index.size()-1, node_labels[node].begin(), node_labels[node].end(), l, r);
		//int block = lead_rs(node + 1) - 1;
		int block = node_blocks[node];

		for (size_type i = l; i <= r; i++) {
			// locate edge
			int occ = index[i];
			int occedge = dels_rs(occ);
			int occedgeindex = occ - ((occedge == 0) ? 0 : dels_ss(occedge) + 1);
			int slength = node_labels[ordered_edges[occedge].first].size();

			// locate specific node in the edge
			int occnode, occnodeindex, occmsaindex;
			if (occedgeindex < slength) { // source node
						      //occnode = std::distance(ordered_node_ids.begin(), std::find(ordered_node_ids.begin(), ordered_node_ids.end(), ordered_edges[occedge].first));
				occnode = ordered_edges[occedge].first;
				//assert(occnode < ordered_node_ids.size());
				occnodeindex = occedgeindex;
			} else { // target node
				 //node = std::distance(ordered_node_ids.begin(), std::find(ordered_node_ids.begin(), ordered_node_ids.end(), ordered_edges[edge].second));
				occnode = ordered_edges[occedge].second;
				//assert(node < ordered_node_ids.size());
				occnodeindex = occedgeindex - slength;
			}

			//int occblock = lead_rs(occnode + 1) - 1;
			int occblock = node_blocks[occnode];

//#ifdef EFG_HPP_DEBUG
//			std::cerr << " edge : " << occedge + 1;
//			std::cerr << " edgeindex : " << occedgeindex + 1;
//			std::cerr << " node : " << occnode + 1;
//			std::cerr << " nodeindex : " << occnodeindex + 1;
//			std::cerr << " block : " << occblock + 1;
//			std::cerr << std::endl;
//#endif

			// semi-repeat-free property (sources and sinks are special)
			if (occnodeindex != 0 || block != occblock) {
#ifdef EFG_HPP_DEBUG
				std::cerr << "Invalid occurrence of node " << node_labels[node] << " (block " << block+1 << ") starting from node " << node_labels[occnode] << " (block " << block+1 << ")" << std::endl;
#endif

				return false;
			}
		}
		return true;
	}

	void efg_validate_worker(
			const size_type startnode,
			const long step,
			const std::vector<std::string> &node_labels,
			const std::vector<std::pair<size_type,size_type>> &ordered_edges,
			const std::vector<size_type> &node_blocks,
			const std::vector<bool> &is_source,
			const std::vector<bool> &is_sink,
			const std::string &ignore_chars,
			const sdsl::csa_wt<> &index,
			const sdsl::rank_support_v5<> &dels_rs,
			const sdsl::select_support_mcl<> &dels_ss,
			std::vector<bool> &to_remove)
	{
		for (size_type i = startnode; i < node_labels.size(); i+= step) {
			if (!efg_validate_node(i, node_labels, ordered_edges, node_blocks, is_source, is_sink, ignore_chars, index, dels_rs, dels_ss)) {
				//to_remove[node_blocks[i]] = true;
				//TODO DEBUG
				assert(node_blocks[i] > 0 and node_blocks[i] < to_remove.size());
				to_remove[node_blocks[i]-1] = true;
			}
		}
	}

	bool efg_validate(
		const std::vector <size_type> &block_indices,
		const std::vector <std::string> &node_labels,
		const std::vector <size_type> &node_blocks,
		const adjacency_list &edges,
		const std::string &ignore_chars,
		const long threads,
		std::vector<bool> &to_remove)
	{
		std::vector<bool> is_source;
		std::vector<bool> is_sink;
		sdsl::bit_vector leaders;
		sdsl::rank_support_v5<> lead_rs;
		sdsl::bit_vector delimiters;
		sdsl::rank_support_v5<> dels_rs;
		sdsl::select_support_mcl<> dels_ss;
		sdsl::csa_wt<> index;

		// build edge index
		std::ostringstream edge_concat;
		size_type cumulative_label_length = 0, cumulative_edge_length = 0, cumulative_walk_length = 0;

		// edge index is concatenation of edge labels with separator chars
		std::vector<std::pair<size_type,size_type>> ordered_edges;
		for (size_type i = 0; i < node_labels.size(); i++) {
			for (size_type j : edges[i]) {
				cumulative_edge_length += node_labels[i].size() + node_labels[j].size();
				ordered_edges.push_back({ i, j });
			}
		}

		delimiters = sdsl::bit_vector(cumulative_edge_length + ordered_edges.size(), 0);

		unsigned long int d = 0;
		for (size_type i = 0; i < node_labels.size(); i++) {
			for (size_type j : edges[i]) {
				edge_concat << node_labels[i] << node_labels[j] << '#';
				d += node_labels[i].size() + node_labels[j].size();
				delimiters[d++] = 1;
			}
		}

//#ifdef EFG_HPP_DEBUG
//		std::cerr << "bitvector delimiters is " << delimiters << std::endl;
//		std::cerr << "edge index is           " << edge_concat.str() << std::endl;
//#endif

		sdsl::construct_im(index, edge_concat.str(), 1);

		// preprocess delimiters
		sdsl::util::init_support(dels_rs, &delimiters);
		sdsl::util::init_support(dels_ss, &delimiters);

		is_source = std::vector<bool>(node_labels.size(), true);
		is_sink = std::vector<bool>(node_labels.size(), true);

		// TODO: parallelize this?
		for (size_type i = 0; i < node_labels.size(); i++) {
			for (size_type j : edges[i]) {
				is_sink[i] = false;
				is_source[j] = false;
			}
		}

		bool result = true;
		if (threads == -1) {
			for (size_type i = 0; i < node_labels.size(); i++) {
				if (!efg_validate_node(i, node_labels, ordered_edges, node_blocks, is_source, is_sink, ignore_chars, index, dels_rs, dels_ss)) {
					//to_remove[node_blocks[i]] = true;
					//TODO DEBUG
					if (node_blocks[i] > 0)
						to_remove[node_blocks[i]-1] = true;
					result = false;
				}
			}
		} else {
			std::vector<std::thread> t;
			for (int i = 0; i < threads; i++) {
				t.push_back(std::thread(efg_validate_worker, i, threads, std::ref(node_labels), std::ref(ordered_edges), std::ref(node_blocks), std::ref(is_source), std::ref(is_sink), std::ref(ignore_chars), std::ref(index), std::ref(dels_rs), std::ref(dels_ss), std::ref(to_remove)));
			}
			for (auto &tt : t)
			{
				tt.join();
			}
			for (size_type i = 0; i < to_remove.size(); i++) {
				if (to_remove[i]) {
					result = false;
					break;
				}
			}
		}

		return result;
	}

}



int main(int argc, char **argv)
{ 
	// TODO: add verbose or log options

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
	const bool elastic(args_info.elastic_flag);
	const bool output_paths(args_info.output_paths_flag);
	const long threads(args_info.threads_arg);
	const bool disable_efg_tricks(args_info.disable_elastic_tricks_flag);

	if (!elastic && output_paths)
	{
		std::cerr << "Output of original sequences as paths without option --elastic is not implemented!\n";
		return EXIT_FAILURE;
	}
	const bool output_gfa_format(args_info.gfa_flag);
	if ((!elastic && output_gfa_format) || (elastic && !output_gfa_format))
	{
		std::cerr << "--elastic and --gfa options are currently only supported when both are used!\n";
		return EXIT_FAILURE;
	}

	if (args_info.heuristic_subset_arg < -1 || args_info.heuristic_subset_arg == 0) {
		std::cerr << "wrong value for --heuristic-subset!\n";
		return EXIT_FAILURE;
	}

	std::string ignorechars;
	if (args_info.ignore_chars_arg != NULL)
		ignorechars = std::string(args_info.ignore_chars_arg);
	else
		ignorechars = "";

	auto start = chrono::high_resolution_clock::now();

	size_type m, n;
	std::vector<std::string> realMSA;
	std::vector<std::string> identifiers;
	if (args_info.heuristic_subset_arg == -1)
	{
		read_input(args_info.input_arg, gap_limit, elastic, realMSA, output_paths, identifiers);

		if (realMSA.empty())
		{
			std::cerr << "Unable to read sequences from the input\n.";
			return EXIT_FAILURE;
		}

		std::cerr << "Input MSA[1.." << realMSA.size() << ",1.." << realMSA[0].size() << "]" << std::endl;   
		m = realMSA.size();
		n = realMSA[0].size();
	} else
	{

		parse_input(args_info.input_arg, m, n, output_paths, identifiers);

		std::cerr << "Input MSA[1.." << m << ",1.." << n << "]" << std::endl;   
	}

	// compute compressed suffix tree of the MSA (or load it if created in a previous run)
	cst_type cst;
	std::vector<std::string> miniMSA;
	size_type miniMSArow = 0;
	if (args_info.heuristic_subset_arg != -1) {
		//miniMSA = std::vector<std::string>(realMSA.begin() + miniMSArow, realMSA.begin() + std::min((long)realMSA.size(), (long)miniMSArow + args_info.heuristic_subset_arg));
		//if (!load_cst(args_info.input_arg, miniMSA, cst, gap_limit))
		//	return EXIT_FAILURE;
		//miniMSArow += args_info.heuristic_subset_arg;
	} else {
		if (!load_cst(args_info.input_arg, realMSA, cst, gap_limit))
			return EXIT_FAILURE;
		std::cerr << "MSA index construction complete, index requires " << sdsl::size_in_mega_bytes(cst) << " MiB." << std::endl;
	}

	std::vector <std::string> node_labels;
	std::vector <size_type> node_blocks; // index of the corresponding block
	std::vector <size_type> block_indices; // starting index of each block
	adjacency_list edges;
	std::vector<std::vector<size_type>> paths;
	std::vector<size_type> f(n, 0); // TODO initialize here
	int status = EXIT_SUCCESS;
	if (elastic) { // semi-repeat-free efg
		if (args_info.heuristic_subset_arg == -1) {
			if (threads == -1)
				segment_elastic_minmaxlength((args_info.heuristic_subset_arg != -1) ? miniMSA : realMSA, cst, ignorechars, block_indices, disable_efg_tricks, f);
			else if (threads > 0)
				segment_elastic_minmaxlength_multithread((args_info.heuristic_subset_arg != -1) ? miniMSA : realMSA, cst, ignorechars, block_indices, threads, disable_efg_tricks, f);
			else {
				std::cerr << "Invalid number of threads." << std::endl;
				return EXIT_FAILURE;
			}
		} else {
			std::cerr << "Starting I/O thread to compute the MSA transpose..." << std::endl;
			std::thread transpose_thread(transpose_msa_worker, args_info.input_arg, m, n, 100); // TODO parameterize buffer
			if (threads == -1) {
				while (miniMSArow < m) {
					miniMSA.clear();
					int _dummy;
					miniMSA = load_rows(args_info.input_arg, args_info.heuristic_subset_arg, _dummy);
					if (!load_cst(args_info.input_arg, miniMSA, cst, gap_limit))
						return EXIT_FAILURE;
					std::cerr << "MSA index construction complete, index requires " << sdsl::size_in_mega_bytes(cst) << " MiB." << std::endl;

					segment_elastic_minmaxlength((args_info.heuristic_subset_arg != -1) ? miniMSA : realMSA, cst, ignorechars, block_indices, disable_efg_tricks, f);
					miniMSArow += args_info.heuristic_subset_arg;
				}
			} else {
				std::vector<std::thread> t;
				for (int i = 0; i < threads; i++) {
					std::cerr << "calling segment_elastic_minmaxlength_worker\n";
					t.push_back(std::thread(segment_elastic_minmaxlength_worker, std::ref(args_info.input_arg), std::ref(ignorechars), std::ref(block_indices), std::ref(disable_efg_tricks), std::ref(f), args_info.heuristic_subset_arg, m, n));
				}

				for (auto &tt : t)
				{
					tt.join();
				}

				segment_elastic_minmaxlength(block_indices, disable_efg_tricks, f, n);
			}
			std::cerr << "Waiting for transpose thread to finish..." << std::flush;
			transpose_thread.join();
			std::cerr << "done." << std::endl;
			fclose(offload_rows_fp);
		}
		f.clear();
		sdsl::util::clear(cst);
	} else if (gap_limit == 1) { // no gaps
		status = segment(realMSA, cst, node_labels, edges);
	} else {
		status = segment2elasticValid(realMSA, cst, node_labels, edges);
	}

	if (status == EXIT_FAILURE)
		return EXIT_FAILURE;

	if (!elastic)
	{
		if (!output_gfa_format) {
			std::cerr << "Writing the index to disk…\n";
			make_efg(block_indices, realMSA, node_labels, node_blocks, edges, output_paths, paths); // TODO: lower RAM usage for huge graphs
			realMSA.clear();
			make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);
		} else {
			std::cerr << "Writing the xGFA to disk…\n";
			make_efg(block_indices, realMSA, node_labels, node_blocks, edges, output_paths, paths); // TODO: lower RAM usage for huge graphs
			realMSA.clear();
			make_gfa(m, n, identifiers, node_labels, edges, node_blocks, block_indices, output_paths, paths, args_info.output_arg, args_info.memory_chart_output_arg);
		}
	}
	if (elastic)
	{
		if (!output_gfa_format) {
			std::cerr << "Writing the index to disk…\n";
			make_efg(block_indices, realMSA, node_labels, node_blocks, edges, output_paths, paths); // TODO: lower RAM usage for huge graphs
			realMSA.clear();
			make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);
		} else {
			if (args_info.heuristic_subset_arg != -1) {
				bool done = false;
				int iterations = 0;
				std::vector<bool> to_remove(block_indices.size(), false);
				while (!done) {
					iterations += 1;
					make_efg_external(args_info.input_arg, m, n, block_indices, realMSA, node_labels, node_blocks, edges, output_paths, paths);
					done = efg_validate(block_indices, node_labels, node_blocks, edges, ignorechars, threads, to_remove);

					std::cerr << "There are ";
					size_type invalidblocks = 0;
					for (auto a : to_remove)
						invalidblocks += a;
					std::cerr << invalidblocks << " blocks to remove" << std::endl;
					/*std::cerr << "the blocks are ";
					for (size_type i = 0; i < block_indices.size(); i++) {
						if (to_remove[i])
							std::cerr << block_indices[i] << " ";
					}
					std::cerr << std::endl;*/

					std::vector <size_type> new_block_indices;
					for (size_type i = 0; i < block_indices.size(); i++) {
						if (!to_remove[i]) {
							new_block_indices.push_back(block_indices[i]);
						} else {
							to_remove[i] = false;
						}
					}
					std::swap(block_indices, new_block_indices);
				}
				std::cerr << "Graph fixed in " << iterations-1 << "iterations…\n";
				std::cerr << "Writing the xGFA to disk…\n";
				output_efg_external(block_indices, args_info.input_arg, m, n, output_paths, identifiers, args_info.output_arg);
			} else {
				std::cerr << "Writing the xGFA to disk…\n";
				output_efg(block_indices, realMSA, output_paths, identifiers, args_info.output_arg);
			}
		}
	}

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::seconds>(end - start);

	if (args_info.graphviz_output_given)
	{
		std::cerr << "Writing the Graphviz file…\n";
		output_graphviz(node_labels, edges, args_info.graphviz_output_arg);
	}

	std::cerr << "Time taken: "
		<< duration.count() << " seconds" << std::endl;

	return EXIT_SUCCESS;
}
