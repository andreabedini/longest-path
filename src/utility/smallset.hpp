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

#ifndef SMALLSET_HPP
#define SMALLSET_HPP

#include <algorithm>
#include <iostream>
#include <vector>

template<class T>
class smallset
{
  std::vector<T> impl_;
public:
  typedef typename std::vector<T>::value_type       value_type;
  typedef typename std::vector<T>::iterator         iterator;
  typedef typename std::vector<T>::const_iterator   const_iterator;
  typedef typename std::vector<T>::reference        reference;
  typedef typename std::vector<T>::const_reference  const_reference;
  typedef typename std::vector<T>::difference_type  difference_type;
  typedef typename std::vector<T>::size_type        size_type;

  smallset() : impl_() {};
  smallset(const smallset& other) : impl_(other.impl_) {}

  iterator begin() { return impl_.begin(); }
  iterator end()   { return impl_.end(); }
  const_iterator begin() const { return impl_.begin(); }
  const_iterator end()   const { return impl_.end(); }

  size_type size() const
  { return impl_.size(); }

  size_type max_size() const
  { return impl_.max_size(); }
  
  bool empty() const
  { return impl_.empty(); }

  void swap(smallset& o)
  { impl_.swap(o.impl_); }

  iterator insert(const T& x) {
    iterator i = std::lower_bound(begin(), end(), x);
    if (i == impl_.end() or x < *i)
      i = impl_.insert(i, x);
    return i;
  }

  void remove(const T& x) {
    iterator i = std::find(begin(), end(), x);
    if (i != impl_.end())
      impl_.erase(i);
  }

  size_type index(const T& x) const {
    return std::distance(begin(), std::find(begin(), impl_.end(), x));
  }

  bool has(const T& x) const {
    return std::binary_search(begin(), end(), x);
  }

  T at(size_type i) const {
    return impl_[i];
  }
};

template<class T>
std::ostream& operator<<(std::ostream& o, const smallset<T>& s) {
    o << "( ";
    const std::size_t n = s.size();
    for (std::size_t i = 0; i < n; ++i)
	o << (int)s.at(i) << " ";
    return o << ")";
}

#endif
