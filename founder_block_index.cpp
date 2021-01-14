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

#include "founder_block_index.hpp"


namespace founder_block_graph {

    auto founder_block_index::serialize(std::ostream &os, sdsl::structure_tree_node* v, std::string name) const -> size_type
    {
        size_type written_bytes(0);
        sdsl::structure_tree_node *child(sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this)));
        written_bytes += csa.serialize(os, child, "csa");
        written_bytes += b_positions.serialize(os, child, "b_positions");
        written_bytes += e_positions.serialize(os, child, "e_positions");
        written_bytes += b_rank1_support.serialize(os, child, "b_rank1_support");
        written_bytes += e_rank1_support.serialize(os, child, "e_rank1_support");
        written_bytes += b_select1_support.serialize(os, child, "b_select1_support");
        written_bytes += e_select1_support.serialize(os, child, "e_select1_support");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }
    
    
    void founder_block_index::load(std::istream& in)
    {
        csa.load(in);
        b_positions.load(in);
        e_positions.load(in);
        b_rank1_support.load(in, &b_positions);
        e_rank1_support.load(in, &e_positions);
        b_select1_support.load(in, &b_positions);
        e_select1_support.load(in, &e_positions);
    }
}
