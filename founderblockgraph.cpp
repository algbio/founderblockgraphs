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
#include "cmdline.h"
#include "founder_block_index.hpp"

/*
  Reads MSA in fasta format as given in --input
  Finds a segment repeat-free segmentation
  Converts the segmentation into a founder block graph
*/

/* To use, install sdsl-lite: https://github.com/simongog/sdsl-lite,
   modify the Makefile, and run make */

/* Copyright (C) 2020 Veli Mäkinen under GNU General Public License v3.0 */


namespace chrono = std::chrono;
namespace fbg = founder_block_graph;

typedef sdsl::cst_sada<> cst_type;
typedef cst_type::node_type node_t;
typedef cst_type::size_type size_type;

typedef std::unordered_set<size_type> edge_set;
typedef std::vector<edge_set> adjacency_list;


class temporary_file
{
protected:
    std::string m_name;
    int         m_fd{-1};
    
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


// For debugging purposes naive string matching
size_type count(std::string const &P,std::string const &T) {
    size_type result = 0;
    for (size_type i=0; i<T.size()-P.size(); i++) 
        if (P==T.substr(i,P.size())) 
            result++;
    return result;  
}

std::string remove_gaps(std:: string const &S) {
    std::string result = "";
    for (size_type i=0; i<S.size(); i++) 
        if (S[i]!='-')
            result += S[i];
    return result;
}

bool check_gaps(std::string const &sequence, std::size_t gap_limit)
{
    if (0 == gap_limit)
        return true;
    
    std::size_t gaprun(0);
    std::size_t maxgaprun(0);
    for (auto const c : sequence)
    {
        if ('-' == c || 'N' == c)
            ++gaprun;
        else
        {
            maxgaprun = std::max(gaprun, maxgaprun);
            gaprun = 0;
        }
    }
    
    maxgaprun = std::max(gaprun, maxgaprun);
    return (maxgaprun < gap_limit);
}


bool check_sequence_length(std::string const &identifier, std::string const &seq, std::size_t expected_length)
{
    if (seq.size() == expected_length)
        return true;
    
    std::cerr << "WARNING: length of the sequence “" << identifier.substr(1) << "” does not match that of the first sequence; skipping. (" << expected_length << " vs. " << seq.size() << ")\n";
    return false;
}


void read_input(char const *input_path, std::size_t gap_limit, std::vector<std::string> &MSA) {
    std::string line, identifier, entry;

    // Reading input fasta
    std::fstream fs;
    fs.open(input_path, std::fstream::in);
    std::getline(fs, identifier); // assuming header first
    
    std::size_t expected_length(0);
    bool is_first(true);
    while (std::getline(fs, line)) {
        if (line[0] == '>') { // header
            if (is_first)
            {
                expected_length = entry.size();
                is_first = false;
            }
            
            if (check_sequence_length(identifier, entry, expected_length) &&
                check_gaps(entry, gap_limit))
            {
                MSA.push_back(entry);
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
        check_gaps(entry, gap_limit))
    {
        MSA.push_back(entry);
    }
}


bool load_cst(char const *input_path, std::vector<std::string> const &MSA, cst_type &cst, size_t gap_limit) {
    
    std::string const index_suffix(".cst");
    std::string const plain_suffix(".plain");
    std::string index_file = std::string(input_path) + plain_suffix + std::to_string(gap_limit) + index_suffix;
    
    // Constructing compressed suffix tree for C
    if (!sdsl::load_from_file(cst, index_file))
    {
        std::cout << "No index "<<index_file<< " located. Building index now." << std::endl;
        
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
            std::cout << "ERROR: File " << input_path << ".plain" << " does not exist. Exit." << std::endl;
            return false;
        }
        sdsl::construct(cst, std::string(input_path)+plain_suffix, 1); // generate index
        sdsl::store_to_file(cst, index_file); // save it
    }
    
    return true;
}


void segment(
    std::vector<std::string> const &MSA,
    cst_type const &cst,
    std::vector <std::string> &out_labels,
    adjacency_list &out_edges
) {

    size_type n = MSA[0].size();
    size_type N = cst.csa.size();

    /* v[j] is such that MSA[1..m][v[j]..j] is a repeat-free segment */
    /* The following computes it naively. 
       This is a preprocessing part of the main algorithm */

    std::vector<size_type> v(n, 0);
    std::unordered_map<size_type,size_type> bwt2row;
    /* bwt2row[(sp,ep)]=i maps BWT[sp..ep] to the row of MSA */

    size_type m = MSA.size();
    std::vector<size_type> sp(m, 0);
    std::vector<size_type> ep(m, cst.csa.size()-1);
    std::vector<size_type> lcp(m, 0);
    std::vector<bool> nested(m,false);
    /* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
    /* lcp[i] = j-jp+1 minus the number of gap symbols in that range*/
    /* nested[i] tells whether current interval is nested in another */

    size_type jp = n; // maintains the left-boundary 
    size_type sum;
    for (size_type j=n-1; j+1>0; j--) {
        v[j] = j+1; // no valid range found
        if (j<n-1) {  
            /* contracting MSA[i][j] from right for all i */
            bwt2row.clear();
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
                bwt2row[sp[i]*N+ep[i]] = i;
            }
        }   
        while (1) {
            bool nestedness = false; 
            // checking nestedness
            // TODO: this takes m^2 time, could be done faster
            for (const auto& pair1 : bwt2row) {
                /* (key=(sa[i],ep[i]), value=i) = (pair.first, pair.second) */              
                nested[pair1.second] = false;
                for (const auto& pair2 : bwt2row)
                    if (pair1!=pair2 && sp[pair1.second]>=sp[pair2.second] && ep[pair1.second]<=ep[pair2.second]) {
                       nested[pair1.second] = true; // extensions to the left may not help  
                       nestedness = true;     
                    }   
            }   
            sum = 0;
            for (const auto& pair : bwt2row)
                /* (key=(sa[i],ep[i]), value=i) = (pair.first, pair.second) */
                if (!nested[pair.second])
                    sum += ep[pair.second]-sp[pair.second]+1;
            if (!nestedness && sum == m) { // no non-aligned matches
                v[j]=jp;
                break;
            }    
            if (nestedness && sum == m) // no other non-aligned matches than nested
                break;
            if (jp==0) // cannot go more to the left
                break;        
            jp--; 
            /* Now looking at MSA[1..m][jp..j] */
            bwt2row.clear();
            for (size_type i=0; i<m; i++) {
                // Keeping old range on gap symbols
                assert(jp < MSA[i].size());
                if (MSA[i][jp]!='-') {
                    sdsl::backward_search(cst.csa,sp[i],ep[i],MSA[i][jp],sp[i],ep[i]);
                    lcp[i]++;
                }
                bwt2row[sp[i]*N+ep[i]]=i;
            }
        } 
    }  

    /* Main algorithm begins */
    /* s[j] is the score of the optimal valid segmentation of MSA[1..m][1..j] */
    /* s[n]=min_{S is valid segmentation} max_{[a..b] \in S} b-a+1 */ 
    std::vector<size_type> s(n, n);

    /* prev[j] is the pointer to the end of previous segment in optimal segmentation */
    std::vector<size_type> prev(n, n);
    // Computation of s[j]'s and prev[j]'s
    for (size_type j=0; j<n; j++) {
        if (v[j]>j) 
            continue; // no valid range
        s[j] = j+1; // handles case jp=0
        prev[j] = j+1; // handles case jp=0
        for (size_type jp=v[j]; jp>0; jp--) {
            if (s[j]>std::max(s[jp-1],j-jp+1)) {
                s[j] = std::max(s[jp-1],j-jp+1);
                prev[j] = jp-1;       
            }
            if (s[j]==j-jp+1) 
                break;
        }
    }
    // outputing optimal score 
    std::cout << "Optimal score: " << s[n-1] << std::endl;

    std::list<size_type> boundariestemp;
    size_type j = n-1;
    boundariestemp.push_front(j);
    while (prev[j]<j) {
        boundariestemp.push_front(prev[j]);
        j = prev[j];         
    }
    std::vector<size_type> boundaries;
    for (const auto& j : boundariestemp)
        boundaries.push_back(j);
    std::cout << "Number of segments: " << boundaries.size() << std::endl;
    //std::cout << "List of segment boundaries: " << std::endl;
    //for (const auto& boundary : boundaries)
    //    std::cout << boundary << std::endl;

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

    std::cout << "#nodes=" << nodecount << std::endl;
    std::cout << "total length of node labels=" << totallength << std::endl; 

    size_type nfounders=0;
    for (size_type j=0; j<boundaries.size(); j++)
        if (blocks[j].size()>nfounders)
            nfounders = blocks[j].size();
    std::cout << "#founders=" << nfounders << std::endl;
/***
    typedef std::unordered_map<size_type, size_type> edge_map;
    typedef std::vector<edge_map> edge_map_vector;
    edge_map_vector edges(nodecount);
    previndex = 0;
    for (size_type k=0; k<boundaries.size()-1; k++) {
        for (size_type i=0; i<m; i++) {
            ellv = remove_gaps(MSA[i].substr(previndex,boundaries[k]-previndex+1));
            ellw = remove_gaps(MSA[i].substr(boundaries[k]+1,boundaries[k+1]-boundaries[k]));  
            edges[str2id[ellv]][str2id[ellw]] = 1;
        }  
*/
    
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
    std::cout << "#edges=" << edgecount << std::endl;
    
    using std::swap;
    swap(out_labels, labels);
    swap(out_edges, edges);
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
    for (std::size_t i(0), count(edges.size()); i < count; ++i)
    {
        auto const &src_label(node_labels[i]);
        
        // Make sure that the destination nodes are sorted.
        auto const &dst_nodes(edges[i]);
        std::vector <size_type> sorted_dst_nodes(dst_nodes.begin(), dst_nodes.end());
        std::sort(sorted_dst_nodes.begin(), sorted_dst_nodes.end());
        
        for (auto const dst_node : sorted_dst_nodes)
        {
            auto const &dst_label(node_labels[dst_node]);
            temp_os << src_label << dst_label << fbg::g_separator_character;
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
    for (auto const &label : node_labels)
    {
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
        assert(match_count);
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


int main(int argc, char **argv) { 

    gengetopt_args_info args_info;
    if (0 != cmdline_parser(argc, argv, &args_info))
        return EXIT_FAILURE;
    
    std::ios_base::sync_with_stdio(false);    // Don't use C style IO after calling cmdline_parser.
    
    auto const gap_limit(args_info.gap_limit_arg);
    if (gap_limit <= 0)
    {
        std::cerr << "Gap limit needs to be non-negative.\n";
        return EXIT_FAILURE;
    }
    
    auto start = chrono::high_resolution_clock::now();

    std::vector<std::string> MSA;
    read_input(args_info.input_arg, gap_limit, MSA);
    std::cout << "Input MSA[1.." << MSA.size() << ",1.." << MSA[0].size() << "]" << std::endl;
    
    cst_type cst;
    if (!load_cst(args_info.input_arg, MSA, cst, gap_limit))
        return EXIT_FAILURE;
    
    std::cout << "Index construction complete, index requires " << sdsl::size_in_mega_bytes(cst) << " MiB." << std::endl;

    std::vector <std::string> node_labels;
    adjacency_list edges;
    segment(MSA, cst, node_labels, edges);
    
    std::cout << "Writing the index to disk…\n";
    make_index(node_labels, edges, args_info.output_arg, args_info.memory_chart_output_arg);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start);
    
    if (args_info.graphviz_output_given)
    {
        std::cout << "Writing the Graphviz file…\n";
        output_graphviz(node_labels, edges, args_info.graphviz_output_arg);
    }
    
    std::cout << "Time taken: "
         << duration.count() << " seconds" << std::endl;
    
    return EXIT_SUCCESS;
}
