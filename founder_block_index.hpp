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

#ifndef FOUNDER_BLOCK_GRAPH_FOUNDER_BLOCK_INDEX_HPP
#define FOUNDER_BLOCK_GRAPH_FOUNDER_BLOCK_INDEX_HPP

#include <istream>
#include <ostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/suffix_array_algorithm.hpp>


namespace founder_block_graph {
    
    typedef sdsl::csa_wt<>  csa_type;
    
    inline constexpr char g_separator_character{'#'};
    
    struct founder_block_index
    {
        typedef csa_type::size_type size_type;
        
        csa_type                    csa;
        sdsl::bit_vector            b_positions;
        sdsl::bit_vector            e_positions;
        sdsl::rank_support_v5<>     b_rank1_support;
        sdsl::rank_support_v5<>     e_rank1_support;
        sdsl::select_support_mcl<>  b_select1_support;
        sdsl::select_support_mcl<>  e_select1_support;
        
        founder_block_index() = default;
        
        founder_block_index(csa_type &&csa_, sdsl::bit_vector &&b_positions_, sdsl::bit_vector &&e_positions_):
            csa(std::move(csa_)),
            b_positions(std::move(b_positions_)),
            e_positions(std::move(e_positions_)),
            b_rank1_support(&b_positions),
            e_rank1_support(&e_positions),
            b_select1_support(&b_positions),
            e_select1_support(&e_positions)
        {
        }
        
        void load(std::istream &in);
        size_type serialize(std::ostream &os, sdsl::structure_tree_node *v = nullptr, std::string name = "") const;
        inline size_type backward_search(char const c, size_type lhs, size_type rhs, size_type &out_lhs, size_type &out_rhs);
        
        template <typename t_reverse_pattern_iterator>
        inline size_type backward_search(t_reverse_pattern_iterator it, t_reverse_pattern_iterator const end);
        
        template <typename t_reverse_pattern_iterator>
        size_type backward_search(t_reverse_pattern_iterator it, t_reverse_pattern_iterator const end, size_type lhs, size_type rhs);
    };
    
    
    auto founder_block_index::backward_search(char const c, size_type lhs, size_type rhs, size_type &out_lhs, size_type &out_rhs) -> size_type
    {
        return sdsl::backward_search(csa, lhs, rhs, c, out_lhs, out_rhs);
    }
    
    
    template <typename t_reverse_pattern_iterator>
    auto founder_block_index::backward_search(t_reverse_pattern_iterator it, t_reverse_pattern_iterator const end) -> size_type
    {
        return backward_search(it, end, 0, csa.size() - 1);
    }


	template <typename t_reverse_pattern_iterator>
	auto founder_block_index::backward_search(t_reverse_pattern_iterator it, t_reverse_pattern_iterator const end, size_type lhs, size_type rhs) -> size_type
	{
		size_type current_count(0);

		while (it != end)
		{
			// Search for the current character.
			std::size_t new_lhs(0);
			std::size_t new_rhs(0);
			current_count = backward_search(*it, lhs, rhs, new_lhs, new_rhs);

			// If not found, try to detect a block pair boundary.
			if (0 == current_count)
			{
				if (!backward_search(g_separator_character, lhs, rhs, new_lhs, new_rhs))
					return 0;

				// Update the lexicographic range.
				// Rank support gives the number of set bits in [0, k).
				auto const r1(b_rank1_support(1 + lhs));
				auto const r2(e_rank1_support(1 + rhs));
				if (! (r1 && r2)) // Sanity check.
					return 0;

				new_lhs = b_select1_support(r1);
				new_rhs = e_select1_support(r2);
				if (! (new_lhs <= lhs && rhs <= new_rhs))
					return 0;

				// Find the current character in the new range.
				current_count = backward_search(*it, new_lhs, new_rhs, lhs, rhs);
			}

			++it;
		}

		return current_count;
	}
}

#endif
