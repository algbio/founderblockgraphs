/**
 * founderblockgraph
 * Copyright (C) 2020 Tuukka Norri
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

#include <climits>
#include <fstream>
#include <iostream>
#include "founder_block_index.hpp"


namespace fbg       = founder_block_graph;
namespace chrono    = std::chrono;


int main(int argc, char **argv)
{
    std::ios_base::sync_with_stdio(false);	// Don't use C style IO.
    
    if (argc < 6)
    {
        std::cerr << "Usage: locate_multiple index_file pattern_file index_id gap_limit repetitions\n";
        return EXIT_FAILURE;
    }
    
    fbg::founder_block_index index;
    
    {
        std::fstream index_is;
        index_is.open(argv[1], std::fstream::in);
        index.load(index_is);
    }
    
    std::string const pattern_file(argv[2]);
    std::string const index_id(argv[3]);
    std::string const gap_limit(argv[4]);
    unsigned long const repetitions(strtoul(argv[5], nullptr, 10));
    if (ULONG_MAX == repetitions)
    {
        std::cerr << "Invalid repetition count: " << strerror(errno) << '\n';
        return EXIT_FAILURE;
    }
    
    std::cout << "INDEX_ID\tGAP_LIMIT\tPATTERN\tLENGTH\tAVERAGE_NS\tDID_FIND\n";
    
    std::string pattern;
    std::size_t pattern_idx(0);
	std::fstream stream;
    stream.exceptions(std::fstream::badbit);
    stream.open(pattern_file, std::ios_base::in);
    while (std::getline(stream, pattern))
    {
        ++pattern_idx;
        //std::cerr << "Pattern " << pattern_idx << "…\n";
        bool is_first(true);
        bool did_find_first(false);
        for (std::size_t i(0); i < repetitions; ++i)
        {
            auto const start(chrono::high_resolution_clock::now());
            founder_block_index::size_type pos(0);
			
            if (is_first)
            {
                did_find_first = (0 == index.backward_search(pattern.begin(), pattern.end(), pos) ? false : true);
                is_first = false;
                
                if (!did_find_first)
                    std::cerr << "WARNING: did not locate “" << pattern << "” in “" << argv[1] << "”, pos = " << pos << ".\n";
            }
            else
            {
                auto const did_find(0 == index.backward_search(pattern.begin(), pattern.end(), pos) ? false : true);
                assert(did_find_first == did_find);
            }
            
            auto const end(chrono::high_resolution_clock::now());
            auto const duration(chrono::duration_cast<chrono::nanoseconds>(end - start));
            auto const current_ms(duration.count());
            
            std::cout << index_id << '\t' << gap_limit << '\t' << pattern_idx << '\t' << pattern.size() << '\t' << current_ms << '\t' << std::uint32_t(did_find_first) << '\n';
        }
    }
    
    return EXIT_SUCCESS;
}
