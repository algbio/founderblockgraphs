#include <iostream>
#include <fstream>      
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <algorithm>
#include <iomanip>
#include <chrono> 

/*
  Reads MSA in fasta format from argv[1]
  Finds a segment repeat-free segmentation
  Converts the segmentation into a founder block graph
*/

/* To use, install sdsl-lite: https://github.com/simongog/sdsl-lite,
   copy this file to its examples subfolder, and run make */

/* Copyright (C) 2020 Veli Mäkinen under GNU General Public License v3.0 */


namespace chrono = std::chrono;

typedef sdsl::cst_sada<>          cst_type;
typedef cst_type::node_type node_t;
typedef cst_type::size_type size_type;


// For debugging purposes naive string matching
size_type count(std::string const &P,std::string const &T) {
    size_type result = 0;
    for (size_type i=0; i<T.size()-P.size(); i++) 
        if (P==T.substr(i,P.size())) 
            result++;
    return result;  
}


bool check_gaps(std::string const &sequence, std::size_t gap_limit)
{
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


void read_input(char const *input_path, std::size_t gap_limit, std::vector<std::string> &MSA) {
    std::string line, entry;

    // Reading input fasta
    std::fstream fs;
    fs.open(input_path, std::fstream::in);
    std::getline(fs,line); // assuming header first
    while (std::getline(fs, line)) {
        if (line[0] == '>') { // header
            if (check_gaps(entry, gap_limit))
                MSA.push_back(entry);
            entry.clear();
        }        
        else
            entry += line;
    }
    fs.close();
    if (check_gaps(entry, gap_limit))
        MSA.push_back(entry); 
}


bool load_cst(char const *input_path, std::vector<std::string> const &MSA, cst_type &cst) {
    
    std::string const index_suffix(".cst");
    std::string const plain_suffix(".plain");
    std::string const index_file(std::string(input_path) + plain_suffix + index_suffix);
    
    // Constructing compressed suffix tree for C
    if (!sdsl::load_from_file(cst, index_file))
    {
        std::cout << "No index "<<index_file<< " located. Building index now." << std::endl;
        
        // Output concatenated inputs to disk. Remove gaps symbols and add separators for indexing.
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


void segment(std::vector<std::string> const &MSA, cst_type const &cst) {

    size_type n = MSA[0].size();

    /* v[j] is such that MSA[1..m][v[j]..j] is a repeat-free segment */
    /* The following computes it naively. 
       This is a preprocessing part of the main algorithm */

    std::vector<size_type> v(n, 0);
    std::unordered_map<size_type, size_type> bwt2row;
    /* bwt2row[j]=i maps j-th smallest suffix to the row of MSA */

    size_type m = MSA.size();
    std::vector<size_type> sp(m, 0);
    std::vector<size_type> ep(m, cst.csa.size()-1);
    std::vector<size_type> lcp(m, 0);
    /* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
    /* lcp[i] = j-jp+1 minus the number of gap symbols in that range*/

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
                    if (cst.depth(u)>=lcp[i]) {
                       sp[i] = cst.lb(u);
                       ep[i] = cst.rb(u);
                    }
                } 
                bwt2row[sp[i]] = i;
            }
        }   
        while (1) {
            sum = 0;
            for (const auto& pair : bwt2row)
                /* (key, value) = (pair.first, pair.second) */
                sum += ep[pair.second]-pair.first+1;
            if (sum == m) { // no non-aligned matches
                v[j]=jp;
                break;
            }    
            if (jp==0)
                break;        
            jp--; 
            /* Now looking at MSA[1..m][jp..j] */
            bwt2row.clear();
            for (size_type i=0; i<m; i++) {
                // Keeping old range on gap symbols
                if (MSA[i][jp]!='-') {
                    sdsl::backward_search(cst.csa,sp[i],ep[i],MSA[i][jp],sp[i],ep[i]);
                    lcp[i]++;
                }
                bwt2row[sp[i]]=i;
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
    for (size_type j=0; j<boundaries.size(); j++) {
        for (size_type i=0; i<m; i++)
            if (!str2id.count(MSA[i].substr(previndex,boundaries[j]-previndex+1))) {       
                blocks[j].push_back(nodecount);
                str2id[MSA[i].substr(previndex,boundaries[j]-previndex+1)] = nodecount++;
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
    
    typedef std::unordered_map<size_type, size_type> edge_map;
    typedef std::vector<edge_map> edge_map_vector;
    edge_map_vector edges(nodecount);
    previndex = 0;
    for (size_type k=0; k<boundaries.size()-1; k++) {
        for (size_type i=0; i<m; i++)
            edges[str2id[MSA[i].substr(previndex,boundaries[k]-previndex+1)]][str2id[MSA[i].substr(previndex,boundaries[k+1]-boundaries[k]+1)]] = 1;
        previndex = boundaries[k]+1;
    }
    size_type edgecount = 0;
    for (size_type i=0; i<nodecount; i++)
        edgecount += edges[i].size();
    std::cout << "#edges=" << edgecount << std::endl;
}


int main(int argc, char **argv) { 
    if (argc <  2) {
        std::cout << "Usage " << argv[0] << " MSA.fasta [GAPLIMIT]" << std::endl;
        std::cout << "    This program constructs a founder block graph of MSA given in FASTA format" << std::endl;
        return 1;
    }
    size_type GAPLIMIT=1; // Filtering out entries with too long gap regions
    if (argc>2)
        GAPLIMIT = atoi(argv[2]);

    auto start = chrono::high_resolution_clock::now();

    std::vector<std::string> MSA;
    read_input(argv[1], GAPLIMIT, MSA);
    std::cout << "Input MSA[1.." << MSA.size() << ",1.." << MSA[0].size() << "]" << std::endl;
    
    cst_type cst;
    if (!load_cst(argv[1], MSA, cst))
        return EXIT_FAILURE;
    
    std::cout << "Index construction complete, index requires " << sdsl::size_in_mega_bytes(cst) << " MiB." << std::endl;

    segment(MSA, cst);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start); 
  
    std::cout << "Time taken: "
         << duration.count() << " seconds" << std::endl;
}
