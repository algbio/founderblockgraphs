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

#include <fstream>
#include <iostream>
#include "founder_block_index.hpp"


namespace fbg = founder_block_graph;


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: locate_partterns index_file\n";
        return EXIT_FAILURE;
    }
    
    fbg::founder_block_index index;
    
    {
        std::fstream index_is;
        index_is.open(argv[1], std::fstream::in);
        index.load(index_is);
    }
    
    while (true)
    {
        std::string pattern;
        
        std::cout << "Pattern? " << std::flush;
        std::cin >> pattern;
        
        if (std::cin.eof())
            return EXIT_SUCCESS;
        
        auto const occurrences(index.backward_search(pattern.rbegin(), pattern.rend()));
        std::cout << occurrences << " occurrences found.\n";
    }
}
