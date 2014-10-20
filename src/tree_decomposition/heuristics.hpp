/*
 *  Created by Andrea Bedini on 24/Nov/2011.
 *  Copyright 2011-2014 Andrea Bedini <andrea@andreabedini.com>
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

#ifndef heuristics_hpp
#define heuristics_hpp

#include <boost/range/algorithm.hpp>

namespace heuristics {
  using namespace boost;

  template<class vertex, class graph>
  void eliminate_vertex(vertex v, graph& g) {
    for (auto a : as_range(adjacent_vertices(v, g))) {
      for (auto b : as_range(adjacent_vertices(v, g))) {
        if (a != b and not edge(a, b, g).second) {
          add_edge(a, b, g);
        }
      }
    }
    clear_vertex(v, g);
    remove_vertex(v, g);
  }

  template<class graph, class vertex>
  unsigned int num_non_adjacent_neighbors(vertex v, graph const& g)
  {
    unsigned int n = 0;
    for (auto u : as_range(adjacent_vertices(v, g))) {
      auto adj = adjacent_vertices(u, g);
      for (auto z : as_range(adjacent_vertices(v, g))) {
        if (boost::find(adj, z) != end(adj))
          n ++;
      }
    }
    return n;
  }

  template<class graph, class outputiterator>
  void greedy_degree_order(graph g, outputiterator out) {
    using vertex = typename graph::vertex_descriptor;
    while (num_vertices(g) > 0) {
      auto v = *min_element(vertices(g), [&](vertex v, vertex u) {
        return degree(v, g) < degree(u, g);
      });
      *out++ = get(vertex_index, g, v);
      eliminate_vertex(v, g);
    }
  }

  //////////////////////////////////////////////////////////////////////

  template<class graph, class outputiterator>
  void greedy_fillin_order(graph g, outputiterator out) {
    using vertex = typename graph::vertex_descriptor;
    while (num_vertices(g) > 0) {
      auto v = *min_element(vertices(g), [&](vertex v, vertex u) {
        return num_non_adjacent_neighbors(v, g) < num_non_adjacent_neighbors(u, g);
      });
      *out++ = get(vertex_index, g, v);
      eliminate_vertex(v, g);
    }
  }

  //////////////////////////////////////////////////////////////////////

  template<class graph, class outputiterator>
  void greedy_local_degree_order(graph g, outputiterator out) {
    using vertex = typename graph::vertex_descriptor;
    auto compare = [&](vertex v, vertex u) {
      return degree(v, g) < degree(u, g);
    };
    vertex current = *min_element(vertices(g), compare);
    while (num_vertices(g) > 0) {
      *out++ = get(vertex_index, g, current);
      vertex next = *min_element(vertices(g), compare);
      eliminate_vertex(current, g);
      current = next;
    }
  }

  //////////////////////////////////////////////////////////////////////

  template<class graph, class outputiterator>
  void greedy_local_fillin_order(graph g, outputiterator out) {
    using vertex = typename graph::vertex_descriptor;
    auto compare =  [&](vertex v, vertex u) {
      return num_non_adjacent_neighbors(v, g) < num_non_adjacent_neighbors(u, g);
    };
    vertex current = *min_element(vertices(g), compare);
    while (num_vertices(g) > 0) {
      *out++ = get(vertex_index, g, current);
      vertex next = *min_element(adjacent_vertices(current, g), compare);
      eliminate_vertex(current, g);
      current = next;
    }
  }
}

#endif
