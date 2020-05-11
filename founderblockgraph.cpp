#include <bits/stdc++.h> 
#include <iostream>
#include <fstream>      
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/suffix_trees.hpp>
#include <algorithm>
#include <iomanip>
#include <chrono> 
using namespace std; 
using namespace sdsl;

/*
  Reads MSA in fasta format from argv[1]
  Finds a segment repeat-free segmentation
  Converts the segmentation into a founder block graph
*/

/* To use, install sdsl-lite: https://github.com/simongog/sdsl-lite,
   copy this file to its examples subfolder, and run make */

/* Copyright (C) 2020 Veli Mäkinen under GNU General Public License v3.0 */


// For debugging purposes naive string matching
size_t count(string P,string T) {
    size_t result = 0;
    for (size_t i=0; i<T.size()-P.size(); i++) 
        if (P==T.substr(i,P.size())) 
            result++;
    return result;  
}

int main(int argc, char **argv) { 
    if (argc <  2) {
        cout << "Usage " << argv[0] << " MSA.fasta [GAPLIMIT]" << endl;
        cout << "    This program constructs a founder block graph of MSA given in FASTA format" << endl;
        return 1;
    }
    size_t GAPLIMIT=1; // Filtering out entries with too long gap regions
    if (argc>2)
        GAPLIMIT = atoi(argv[2]);

    auto start = chrono::high_resolution_clock::now();

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


    // Filtering out entries with too many gaps or N's
    size_t ngaps=0;
    size_t maxgaprun=0;
    size_t gaprun=0;
    for (size_t i=0; i<MSAtemp.size(); i++) {
        for (size_t j=0; j<MSAtemp[i].size(); j++)
            if (MSAtemp[i][j]=='-' || MSAtemp[i][j]=='N') {
               ngaps++;
               gaprun++;
            } 
            else {
               if (gaprun>maxgaprun)
                   maxgaprun = gaprun;
               gaprun = 0; 
            }
        if (gaprun>maxgaprun)
           maxgaprun = gaprun;
        gaprun = 0;
        if (maxgaprun<GAPLIMIT) 
            MSA.push_back(MSAtemp[i]);
        ngaps = 0;
        maxgaprun = 0;  
    }     
    cout << "Input MSA[1.." << MSA.size() << ",1.." << MSA[0].size() << "]" << endl; 
    // Concatenating, removing gap symbols, and adding separators for indexing
    string C="";
    for (size_t i=0; i<MSA.size(); i++) {
        for (size_t j=0; j<MSA[i].size(); j++) 
            if (MSA[i][j]!='-') 
                C += MSA[i][j]; 
        C += '#';
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
       This is a preprocessing part of the main algorithm */

    size_t *v = new size_t[n];
    for (size_t i=0; i<n; i++)
        v[i] = 0;
 
    unordered_map<size_t, size_t> bwt2row;
    /* bwt2row[j]=i maps j-th smallest suffix to the row of MSA */

    size_t m = MSA.size();
    size_t *sp = new size_t[m];
    size_t *ep = new size_t[m];
    size_t *lcp = new size_t[m];
    /* [sp[i]..ep[i]] maintains BWT interval of MSA[i][jp..j] */
    /* lcp[i] = j-jp+1 minus the number of gap symbols in that range*/

    // initializing
    for (size_t i=0; i<m; i++) {
        sp[i] = 0;
        ep[i] = cst.csa.size()-1; 
        lcp[i] = 0;
    }   

    size_t jp = n; // maintains the left-boundary 
    size_t sum;
    for (size_t j=n-1; j+1>0; j--) {
        v[j] = j+1; // no valid range found
        if (j<n-1) {  
            /* contracting MSA[i][j] from right for all i */
            bwt2row.clear();
            for (size_t i=0; i<m; i++) {
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
            for (size_t i=0; i<m; i++) {
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

    list<size_t> boundariestemp;
    size_t j = n-1;
    boundariestemp.push_front(j);
    while (prev[j]<j) {
        boundariestemp.push_front(prev[j]);
        j = prev[j];         
    }
    vector<size_t> boundaries;
    for (const auto& j : boundariestemp)
        boundaries.push_back(j);
    cout << "Number of segments: " << boundaries.size() << endl;
    //cout << "List of segment boundaries: " << endl;
    //for (const auto& boundary : boundaries)
    //    cout << boundary << endl;

    /* Convert segmentation into founder block graph */
    unordered_map<string, size_t> str2id;
    size_t nodecount = 0; 
    size_t previndex = 0;
    vector<size_t> *blocks = new vector<size_t>[boundaries.size()];
    for (size_t j=0; j<boundaries.size(); j++) {
        for (size_t i=0; i<m; i++)
            if (!str2id.count(MSA[i].substr(previndex,boundaries[j]-previndex+1))) {       
                blocks[j].push_back(nodecount);
                str2id[MSA[i].substr(previndex,boundaries[j]-previndex+1)] = nodecount++;
            }
        previndex = boundaries[j]+1;
    }

    string *labels = new string[nodecount];
    for (const auto& pair : str2id) {
        labels[pair.second] = pair.first;
    }
    size_t totallength = 0; 
    for (size_t i=0; i<nodecount; i++)
        totallength += labels[i].size();

    cout << "#nodes=" << nodecount << endl;
    cout << "total length of node labels=" << totallength << endl; 
    
    unordered_map<size_t, size_t> *edges = new unordered_map<size_t, size_t>[nodecount];
    previndex = 0;
    for (size_t k=0; k<boundaries.size()-1; k++) {
        for (size_t i=0; i<m; i++)
            edges[str2id[MSA[i].substr(previndex,boundaries[k]-previndex+1)]][str2id[MSA[i].substr(previndex,boundaries[k+1]-boundaries[k]+1)]] = 1;
        previndex = boundaries[k]+1;
    }
    size_t edgecount = 0;
    for (size_t i=0; i<nodecount; i++)
        edgecount += edges[i].size();
    cout << "#edges=" << edgecount << endl;

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - start); 
  
    cout << "Time taken: "
         << duration.count() << " seconds" << endl;     
}


