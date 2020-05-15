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
        sdsl::select_support_mcl<>  b_select1_support;
        sdsl::select_support_mcl<>  e_select1_support;
        
        founder_block_index() = default;
        
        founder_block_index(csa_type &&csa_, sdsl::bit_vector &&b_positions_, sdsl::bit_vector &&e_positions_):
            csa(std::move(csa_)),
            b_positions(std::move(b_positions_)),
            e_positions(std::move(e_positions_)),
            b_rank1_support(&b_positions),
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
        size_type retval(0);
        size_type current_count(0);
        
        if (it != end)
        {
            // Search for the first character.
            current_count = backward_search(*it, lhs, rhs, lhs, rhs);
            if (0 == current_count)
                return 0;
            
            ++it;
            if (it == end)
                return current_count;
            else
            {
                do
                {
                    // Try to detect a block pair boundary.
                    // To this end, branch with g_separator_character on each step.
                    {
                        size_type new_lhs(0);
                        size_type new_rhs(0);
                        current_count = backward_search(g_separator_character, lhs, rhs, new_lhs, new_rhs);
                        if (current_count)
                        {
                            // Update the lexicographic range.
                            // Rank support gives the number of set bits in [0, k).
                            auto const r(b_rank1_support(1 + lhs));
                            if (r)
                            {
                                new_lhs = b_select1_support(r);
                                new_rhs = e_select1_support(r);
                                if (new_lhs <= lhs && rhs <= new_rhs)
                                    retval += backward_search(it, end, new_lhs, new_rhs);
                            }
                        }
                    }
                    
                    // Continue the search in the current block.
                    current_count = backward_search(*it, lhs, rhs, lhs, rhs);
                    if (0 == current_count)
                        break;
                    ++it;
                } while (it != end);
                retval += current_count;
            }
        }
        
        return retval;
    }
}

#endif
