/*
 *  Created by Andrea Bedini on 27/Feb/2013.
 *  Copyright 2013-2014 Andrea Bedini <andrea@andreabedini.com>
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

#ifndef LONGEST_PATH_HPP
#define LONGEST_PATH_HPP

#include "boost/format.hpp"
#include "boost/math/tools/polynomial.hpp"
#include "boost/range/algorithm/count.hpp"
#include "boost/range/algorithm/find_if.hpp"
#include "boost/range/algorithm/replace.hpp"
#include "boost/unordered/unordered_map.hpp"

#include "boost/optional.hpp"

template<class Weight>
struct longest_path
{
  using weight_type = Weight ;
  using connectivity = std::vector<int8_t>;
  using table_type = boost::unordered_map<connectivity, weight_type>;

  static bool is_endpoint(connectivity const& c, size_t i)
  {
    return c[i] > 0 and boost::count(c, c[i]) == 1;
  }

  static boost::optional<connectivity>
  connect(connectivity c, size_t i, size_t j)
  {
    // don't touch a finished state
    if (is_finished(c))
      return {};

    auto li = c[i], lj = c[j];

    // check we are not hitting bullets
    if (li < 0 or lj < 0)
      return {};

    // they are both empty strands
    if (li == 0 and lj == 0) {
      c[i] = c[j] = *boost::max_element(c) + 1;
      return c;
    }

    // i is not empty, j is empty
    if (li != 0 and lj == 0) {
      c[i] = -1; c[j] = li;
      return c;
    }

    // i is empty, j is not empty
    if (li == 0 and lj != 0) {
      c[i] = lj; c[j] = -1;
      return c;
    }

    // i and j both non empty
    assert(li != 0 and lj != 0);

    // check we are not closing a loop
    if (li == lj)
      return {};

    // check other strands
    if (is_endpoint(c, i) and is_endpoint(c, j)) {
      // check if there are no non-empty strands beside i and j
      if (boost::count_if(c, [=](int x) { return x > 0 and x != li and x != lj; }) == 0) {
        // in this case we have a finished state
        c.clear();
        return c;
      } else {
        // otherwise the state is invalid
        return {};
      }
    }

    boost::replace(c, lj, li);
    c[i] = c[j] = -1;
    return c;
  }

  static connectivity detach(connectivity c, size_t i)
  {
    c[i] = *boost::max_element(c) + 1;
    return c;
  }

  static connectivity canonicalize(connectivity c)
  {
    const auto max = 50;
    int8_t table[max];
    boost::fill(table, 0);

    // starts recounting from 1
    auto k = 1;
    for (auto& x : c) {
      if (x <= 0)
        continue;
      if (table[x] != 0) {
        x = table[x];
      } else {
        table[x] = k;
        x = k;
        ++ k;
      }
    }
    return c;
  }

  static size_t how_many_endpoints(connectivity const& c)
  {
    std::vector<int> table(c.size() + 1, -1);
    size_t result = 0;
    for (size_t j = 0; j < c.size(); ++j) {
      if (c[j] > 0) {
        if (table[c[j]] < 0) {
          table[c[j]] = j;
          result ++;
        } else {
          result --;
        }
      }
    }
    return result;
  }

  static bool is_empty(connectivity const& c)
  {
    return not c.empty() and boost::count_if(c, [](int x) { return x != 0; }) == 0;
  }

  static bool is_finished(connectivity const& c)
  {
    return c.empty();
  }

  template<class F>
  static void decompose(connectivity const& c, F f)
  {
    std::map<size_t, size_t> table;

    auto j = 0;
    for (auto& x : c) {
      if (x <= 0)
        continue;
      if (table[x] != -1)
        f(table[x], j);
      table[x] = j;
      j ++;
    }
  }

  static table_type empty_state(size_t size)
  {
    return table_type{ { connectivity(size), weight_type(1) } };
  }

  table_type
  join_operator(size_t i, size_t j, table_type const& table) const
  {
    table_type new_table{table};
    for (auto const& state : table) {
      auto maybe_newc = connect(state.first, i, j);
      if (maybe_newc and how_many_endpoints(*maybe_newc) <= 2) {
        new_table[canonicalize(*maybe_newc)] += (state.second << 1);
      }
    }
    return new_table;
  }

  static boost::optional<connectivity>
  delete_node(connectivity const& c, size_t i)
  {
    if (is_finished(c))
      return c;

    if (is_endpoint(c, i)) {
        // if there are other strands
        if (boost::count_if(c, [&](int x) { return x > 0 and x != c[i]; }) == 0) {
          // if there are no other strands
          return connectivity();
        } else {
          return {};
        }
    }
    connectivity newc{c};
    newc.erase(newc.begin() + i);
    return canonicalize(newc);
  }

  table_type
  delete_operator(size_t i, table_type const& table) const
  {
    table_type new_table;
    for (auto const& state : table) {
      auto maybe_newc = delete_node(state.first, i);
      if (maybe_newc and how_many_endpoints(*maybe_newc) <= 2)
        new_table[*maybe_newc] += state.second;
    }
    return new_table;
  }

  template<class Mapping>
  table_type
  table_fusion(Mapping A_to_B, table_type const& A_table, table_type const& B_table) const
  {
    table_type new_table;
    for (auto const& stateA : A_table) {
      for (auto const& stateB : B_table) {
        // if state A is finished
        if (is_finished(stateA.first)) {
          if (is_empty(stateB.first)) {
            new_table[connectivity()] += stateA.second * stateB.second;
          }
          continue;
        }

        // if state B is finished
        if (is_finished(stateB.first)) {
          if (is_empty(stateA.first)) {
            new_table[connectivity()] += stateA.second * stateB.second;
          }
          continue;
        }

        // n is the size of the destination (stateB)
        auto const n = stateB.first.size();

        // convert to the new order
        connectivity newa(n);
        for (auto i = 0; i < stateA.first.size(); ++i)
          newa[A_to_B[i]] = stateA.first[i];

        connectivity newc(stateB.first);

        bool valid = true;
        std::map<size_t, size_t> table;
        for (size_t i = 0; i < n; ++i) {
          if (newa[i] > 0) {
            // a strand connected or not
            if (table.find(newa[i]) != table.end()) {
              // we have a link
              auto j = table[newa[i]];
              // go ahead and connect
              if (auto maybe_newc = connect(newc, i, j)) {
                newc = *maybe_newc;
              } else {
                valid = false;
                break;
              }
              table.erase(newa[i]);
            } else {
              table[newa[i]] = i;
            }
          }
          // a bullet
          if (newa[i] == -1) {
            if (newc[i] != 0) {
              valid = false;
              break;
            }
            newc[i] = -1;
          }
        }
        if (not valid) {
          continue;  // to exit the 'for' (to the next pair of states)
        }
        // now the single strands
        for (auto x : table) {
          // discard the state if we have single strands to apply to a finished state
          if (is_finished(newc)) {
            valid = false;
            break;
          }
          auto bi = x.second;
          switch (newc[bi]) {
            case -1:
              valid = false;
              break;
            case 0:
              newc = detach(newc, bi);
              break;
            default:
              // we can do this only if there is nothing else
              if (boost::count(newc, newc[bi]) == 2) {
                // part of a pair
                newc[bi] = -1;
              } else if (boost::count_if(newc, [](int x) { return x > 0; }) == 1) {
                // single strand
                newc.clear(); // mark as finished
              } else {
                // there's other stuff
                valid = false;
                continue;
              }
          }
        }
        if (valid and how_many_endpoints(newc) <= 2) {
          new_table[canonicalize(newc)] += stateA.second * stateB.second;
        }
      }
    }
    return new_table;
  }
};

#endif
