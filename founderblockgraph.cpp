#include <iostream>
#include <fstream>
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
#include <vector>
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

namespace {

	namespace chrono = std::chrono;
	namespace fbg = founder_block_graph;

	typedef sdsl::cst_sct3<> cst_type;
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

		// Opening output stream for DEBUGing
		std::fstream fs2;
		fs2.open((std::string(input_path)+"."+std::to_string(gap_limit) + ".filtered").c_str(), std::fstream::out);

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
					fs2 << remove_gaps(entry) << std::endl;
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
			fs2 << remove_gaps(entry) << std::endl;
		}
		fs2.close();      
	} 


	bool load_cst(char const *input_path, std::vector<std::string> const &MSA, cst_type &cst, size_t gap_limit) {

		std::string const index_suffix(".cst");
		std::string const plain_suffix(".plain");
		std::string index_file = std::string(input_path) + plain_suffix + std::to_string(gap_limit) + index_suffix;

		// Constructing compressed suffix tree for C
		if (!sdsl::load_from_file(cst, index_file))
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

			std::ifstream in(std::string(input_path)+plain_suffix);
			if (!in) {
				std::cerr << "ERROR: File " << input_path << ".plain" << " does not exist. Exit." << std::endl;
				return false;
			}
			sdsl::construct(cst, std::string(input_path)+plain_suffix, 1); // generate index
			sdsl::store_to_file(cst, index_file); // save it
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

	void segment_elastic_minmaxlength(
			std::vector<std::string> const &MSA,
			cst_type const &cst,
			std::string &ignorechars,
			std::vector<size_type> &out_indices
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
		sdsl::bit_vector fullrow(m, true); // mark if leaves[i] is still the initial value

		/* binary coloring of the leaves */
		sdsl::bit_vector color(cst.size(cst.root()), false);

		/* f[x] is minimum index greater or equal to x such that MSA[0..m-1][x..f[x]] is a semi-repeat-free segment */
		std::vector<size_type> f(n, 0);
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
									fi = indexedrows_ss[ii].select(indexedrows_rs[ii].rank(n));
									// fi can be less than x here but it's still correct
								} else {
									fi = indexedrows_ss[ii].select(indexedrows_rs[ii].rank(x) + g);
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
			f[x] = fimax;

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

		// there is always a valid segmentation
		/*if (f[0] == n) {
			std::cerr << "No valid segmentation found!\n";
			return ;
		}*/
		/*std::cerr << "f : ";
		for (auto v : f)
			std::cerr << v << " ";
		std::cerr << std::endl;*/

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
			while (j == std::get<1>(minimal_right_extensions[y])) {
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

				if (should_stop)
					std::abort();
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
			std::vector<std::string> const &MSA,
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
		dst_os << "M\t" << MSA.size() << "\t" << MSA[0].length() << std::endl;

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
			assert(MSA.size() == paths.size());
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
	if (!elastic && output_paths)
	{
		std::cerr << "Output of original sequences as paths without option --elastic is not implemented!\n";
		return EXIT_FAILURE;
	}
	const bool output_gfa_format(args_info.gfa_flag);
	if (!elastic && output_gfa_format)
	{
		std::cerr << "xGFA output without option --elastic not yet implemented!\n";
		return EXIT_FAILURE;
	}

	std::string ignorechars;
	if (args_info.ignore_chars_arg != NULL)
		ignorechars = std::string(args_info.ignore_chars_arg);
	else
		ignorechars = "";

	auto start = chrono::high_resolution_clock::now();

	std::vector<std::string> MSA;
	std::vector<std::string> identifiers;
	read_input(args_info.input_arg, gap_limit, elastic, MSA, output_paths, identifiers);
	if (MSA.empty())
	{
		std::cerr << "Unable to read sequences from the input\n.";
		return EXIT_FAILURE;
	}
	std::cerr << "Input MSA[1.." << MSA.size() << ",1.." << MSA[0].size() << "]" << std::endl;   

	// compute compressed suffix tree of the MSA (or load it if created in a previous run)
	cst_type cst;
	if (!load_cst(args_info.input_arg, MSA, cst, gap_limit))
		return EXIT_FAILURE;

	std::cerr << "MSA index construction complete, index requires " << sdsl::size_in_mega_bytes(cst) << " MiB." << std::endl;

	std::vector <std::string> node_labels;
	std::vector <size_type> node_blocks; // index of the corresponding block
	std::vector <size_type> block_indices; // starting index of each block
	adjacency_list edges;
	std::vector<std::vector<size_type>> paths;
	int status = EXIT_SUCCESS;
	if (elastic) { // semi-repeat-free efg
		segment_elastic_minmaxlength(MSA, cst, ignorechars, block_indices);
		make_efg(block_indices, MSA, node_labels, node_blocks, edges, output_paths, paths);
	} else if (gap_limit == 1) { // no gaps
		status = segment(MSA, cst, node_labels, edges);
	} else {
		status = segment2elasticValid(MSA, cst, node_labels, edges);
	}

	if (status == EXIT_FAILURE)
		return EXIT_FAILURE;

	if (!elastic)
	{
		if (!output_gfa_format) {
			std::cerr << "Writing the index to disk…\n";
			make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);
		} else {
			std::cerr << "Writing the xGFA to disk…\n";
			make_gfa(MSA, identifiers, node_labels, edges, node_blocks, block_indices, output_paths, paths, args_info.output_arg, args_info.memory_chart_output_arg);
		}
	}
	if (elastic)
	{
		if (!output_gfa_format) {
			std::cerr << "Writing the index to disk…\n";
			make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);
		} else {
			std::cerr << "Writing the xGFA to disk…\n";
			make_gfa(MSA, identifiers, node_labels, edges, node_blocks, block_indices, output_paths, paths, args_info.output_arg, args_info.memory_chart_output_arg);
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
