#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>      
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <algorithm>
#include <iomanip>
using namespace std; 
using namespace sdsl;

#define GAPLIMIT 50 // Filtering out entries with too long gap regions

/*
# reads MSA in fasta format from argv[1]
# finds a block repeat-free segmentation
# converts the segmentation into a founder block graph
*/

/* To use, install sdsl-lite: https://github.com/simongog/sdsl-lite,
   copy this file to its examples subfolder, and run make */


// For debugging purposes naive string matching
unsigned int count(string P,string T) {
    unsigned int result = 0;
    for (unsigned int i=0; i<T.size()-P.size(); i++) 
        if (P==T.substr(i,P.size())) 
            result++;
    return result;  
}

int main(int argc, char **argv) { 
    if (argc <  2) {
        cout << "Usage " << argv[0] << " MSA.fasta " << endl;
        cout << "    This program constructs a founder block graph of MSA given in FASTA format" << endl;
        return 1;
    }
    string line,entry;
    vector<string> MSAtemp;
    vector<string> MSA;

    // Reading input fasta
    std::fstream fs;
    fs.open(argv[1], fstream::in);  
    getline(fs,line); // assuming header first
    while (getline(fs,line)) {
        if (line[0]=='>') { // header            
            MSAtemp.push_back(entry);  
            entry = "";     
        }        
        else
            entry += line;
    }
    fs.close(); 
    MSAtemp.push_back(entry); 


    // Filtering out entries with too many gaps
    size_t ngaps=0;
    size_t maxgaprun=0;
    size_t gaprun=0;
    for (size_t i=0; i<MSAtemp.size(); i++) {
        for (size_t j=0; j<MSAtemp[i].size()-1; j++)
            if (MSAtemp[i][j]=='-') {
               ngaps++;
               if (MSAtemp[i][j+1]=='-') 
                  gaprun++; 
            } 
            else {
               if (gaprun>maxgaprun)
                   maxgaprun = gaprun;
               gaprun = 1; 
            }
        if (gaprun>maxgaprun)
           maxgaprun = gaprun;     
        if (maxgaprun<GAPLIMIT) 
            MSA.push_back(MSAtemp[i]);
        ngaps = 0;
        maxgaprun = 0;  
    }     
    cout << "Input MSA[1.." << MSA.size() << ",1.." << MSA[0].size() << "]" << endl; 
    // Concatenating for indexing
    string C="";
    for (size_t i=0; i<MSA.size();i++) {
        C += MSA[i];
    }

    // Outputing concatenation to disk
    string plain_suffix = ".plain"; 
    fs.open(string(argv[1]) + plain_suffix, fstream::out);  
    fs << C;
    fs.close();

    // Constructing compressed suffix tree for C
    string index_suffix = ".cst";
    
    string index_file   = string(argv[1])+plain_suffix+index_suffix;

    typedef cst_sada<>::node_type node_t;
    cst_sada<> cst;
    if (!load_from_file(cst, index_file)) {
        ifstream in(string(argv[1])+plain_suffix);
        if (!in) {
            cout << "ERROR: File " << argv[1] << ".plain" << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index "<<index_file<< " located. Building index now." << endl;
        construct(cst, string(argv[1])+plain_suffix, 1); // generate index
        store_to_file(cst, index_file); // save it
    }
    cout << "Index construction complete, index requires " << size_in_mega_bytes(cst) << " MiB." << endl;


    size_t n = MSA[0].size();

    /* v[j] is such that MSA[1..m][v[j]..j] is a repeat-free segment */
    /* The following computes it naively. 
       This is a preprocessing part of the main algorithm*/

    size_t *v = new size_t[n];
    for (size_t i=0; i<n; i++)
        v[i] = 0;
 
    unordered_map<size_t, size_t> ranges;
    /* ranges[2]=5 tells that BWT interval [2..5] contains matches 
       for subset of rows MSA[1..m][jp..j] */
    size_t m = MSA.size();
    size_t *sp = new size_t[m];
    size_t *ep = new size_t[m];
    /* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
    // initializing
    for (size_t i=0; i<m; i++) {
        sp[i] = 0;
        ep[i] = cst.csa.size()-1; 
    }   

    size_t jp = n; // maintains the left-boundary 
    size_t sum;
    for (size_t j=n-1; j+1>0; j--) {
        v[j] = j+1; // no valid range found
        //cout << "[" << jp << ".." << j << "]" << endl;
        if (j<n-1) {  
            /* contracting MSA[i][j] from right for all i */
            ranges.clear();
            for (size_t i=0; i<m; i++) {
                /* mapping BWT range to suffix tree node,
                   taking parent, mapping back if parent is not too far*/ 
                node_t ll = cst.select_leaf(sp[i]+1); 
                node_t rl = cst.select_leaf(ep[i]+1); 
                node_t w = cst.lca(ll,rl);
                node_t u = cst.parent(w);
                if (cst.depth(u)>=j-jp+1) {
                   sp[i] = cst.lb(u);
                   ep[i] = cst.rb(u);
                }
                ranges[sp[i]] = ep[i];
            }
        }   
        while (1) {
            sum = 0;
            for (const auto& range : ranges)
                /* (key, value) = (range.first, range.second) */
                sum += range.second-range.first+1;
            //if (sum<m and sum>0) {
            //    cout << "this should never happen: sum < m: " << sum << "<" << m << endl;
            //} 
            if (sum == m) { // no non-aligned matches
                v[j]=jp;
                // cout << "v[" << j << "]=" << jp << endl;
                break;
            }    
            if (jp==0)
                break;        
            jp--; 
            /* Now looking at MSA[1..m][jp..j] */
            ranges.clear();
            for (size_t i=0; i<m; i++) {
                sdsl::backward_search(cst.csa,sp[i],ep[i],MSA[i][jp],sp[i],ep[i]);
                ranges[sp[i]]=ep[i];
            }
        } 
    }  

    /* Main algorithm begins */
    /* s[j] is the score of the optimal valid segmentation of MSA[1..m][1..j] */
    /* s[n]=min_{S is valid segmentation} max_{[a..b] \in S} b-a+1 */ 
    size_t *s = new size_t[n];

    /* prev[j] is the pointer to the end of previous segment in optimal segmentation */
    size_t *prev = new size_t[n];
    // Computation of s[j]'s and prev[j]'s
    for (size_t j=0; j<n; j++) {
        s[j] = n; // init
        prev[j] = n; // init
        if (v[j]>j) 
            continue; // no valid range
        s[j] = j+1; // handles case jp=0
        prev[j] = j+1; // handles case jp=0
        for (size_t jp=v[j]; jp>0; jp--) {
            if (s[j]>max(s[jp-1],j-jp+1)) {
                s[j] = max(s[jp-1],j-jp+1);
                prev[j] = jp-1;       
            }
            if (s[j]==j-jp+1) 
                break;
        }
    }
    // outputing optimal score 
    cout << "Optimal score: " << s[n-1] << endl;

    list<size_t> boundaries;
    size_t j = n-1;
    boundaries.push_front(j);
    while (prev[j]<j) {
        boundaries.push_front(prev[j]);
        j = prev[j];         
    }
    cout << "Number of segments: " << boundaries.size() << endl;
    cout << "List of segment boundaries: " << endl;
    for (const auto& boundary : boundaries)
        cout << boundary << endl;
    /* TODO: Convert segmentation into founder block graph */

}


