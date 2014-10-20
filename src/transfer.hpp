/*
 *  Created by Andrea Bedini on 13/Jul/2009.
 *  Copyright 2009-2014 Andrea Bedini <andrea@andreabedini.com>
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */

#ifndef TRANSFER_HPP
#define TRANSFER_HPP

#include "tree_decomposition/tree_decomposition.hpp"
#include <boost/range/algorithm/set_algorithm.hpp>

namespace transfer {
  using tree_decomposition::vertex_list;
  using tree_decomposition::bag_ptr;

  template<class Operators>
  typename Operators::table_type
  recurse(const Operators& op, bag_ptr b)
  {
    // create a new table containing only the empty state
    auto const n = b->vertices.size();
    auto table = op.empty_state(n);

    // iterates over children
    for (auto b_sib : b->children) {
      // recurse
      auto table_sib = recurse(op, b_sib);

      // diffe contains the vertices in b_sib which are not in b (the parent bag)
      std::vector<unsigned int> diffe;
      boost::set_difference(b_sib->vertices, b->vertices, std::back_inserter(diffe));

      // delete each vertex not present in the parent bag

      // we need to make a copy first because
      // 1) we need to keep the indices consistent while removing vertices
      // 2) we don't want to destroy the tree decomposition
      vertex_list b_sib_left_over(b_sib->vertices);
      for (auto v : diffe) {
        table_sib = op.delete_operator(b_sib_left_over.index(v), table_sib);
        b_sib_left_over.remove(v);
      }

      // create b_sib to b bag mapping
      auto const A_size = b_sib_left_over.size();
      std::vector<unsigned int> A_to_B(A_size);
      for (unsigned int i = 0; i < A_size; ++i)
        A_to_B[i] = b->vertices.index(b_sib_left_over.at(i));

      table = op.table_fusion(A_to_B, table_sib, table);
    }

    // apply the join operator for each edge in the bag
    for (auto e : b->edges) {
      table = op.join_operator(b->vertices.index(e.first),
        b->vertices.index(e.second), table);
    }
    return table;
  }

  template<class Operators>
  typename Operators::weight_type
  transfer(const Operators& op, bag_ptr b)
  {
    auto table = recurse(op, b);

    // we need to make a copy first because
    // 1) we need to keep the indices consistent while removing vertices
    // 2) we don't want to destroy the tree decomposition
    vertex_list v_to_remove(b->vertices);
    for (auto v : b->vertices) {
      table = op.delete_operator(v_to_remove.index(v), table);
      v_to_remove.remove(v);
    }
    assert(table.size() == 1);
    return table.begin()->second;
  }
}

#endif
