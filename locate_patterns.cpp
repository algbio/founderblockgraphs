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
#include "locate_patterns_cmdline.h"


namespace fbg = founder_block_graph;


int main(int argc, char **argv)
{
    gengetopt_args_info args_info;
    if (0 != cmdline_parser(argc, argv, &args_info))
        return EXIT_FAILURE;
    
    std::ios_base::sync_with_stdio(false);    // Don't use C style IO after calling cmdline_parser.
    
    fbg::founder_block_index index;
    
    {
        std::fstream index_is;
        index_is.open(args_info.index_arg, std::fstream::in);
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

		if (args_info.error_on_not_found_flag && 0 == occurrences)
		{
			std::cerr << "ERROR: pattern not found.\n";
			return EXIT_FAILURE;
		}
    }
}
