/*
 *  Created by Andrea Bedini on 21/Nov/2008.
 *  Copyright 2008-2014 Andrea Bedini <andrea@andreabedini.com>
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

#ifndef POLYNOMIAL_HPP_
#define POLYNOMIAL_HPP_

#include <boost/operators.hpp>

#include <iosfwd>
#include <vector>

template<typename T>
class polynomial
  : boost::ring_operators1< polynomial<T>
  , boost::ring_operators2< polynomial<T>, T
  , boost::equality_comparable< polynomial<T>
  , boost::left_shiftable<polynomial<T>, unsigned int
  > > > >
{
  /// Coefficients { c_0; ...; c_n } :
  std::vector<T> impl_;

  /// remove leading zero coefficients
  void normalize()
  {
    std::size_t i = order();
    // do not remove the constant term, even if it's zero
    while (impl_[i] == 0 and i > 1)
      i--;
    order(i);
  }

public:
  // default constructor
  explicit polynomial(const T& a = T(0))
    : impl_{a}
  {
  }

  polynomial(std::initializer_list<T> list)
    : impl_(list)
  {
  }

  // copy constructor
  polynomial(polynomial<T> const& rhs)
    : impl_(rhs.impl_)
  {
  }

  // move constructor
  polynomial(polynomial<T>&& rhs)
    : impl_(std::move(rhs.impl_))
  {
  }
  
  // assignemnt
  polynomial<T>& operator=(polynomial<T> const& rhs)
  {
    impl_ = rhs.impl_;
    return *this;
  }

  polynomial<T>& operator=(polynomial<T>&& rhs)
  {
    impl_ = std::move(rhs.impl_);
    return *this;
  }

  // conversions
  template<typename T2>
  friend class polynomial;

  template<class T2>
  explicit polynomial(T2 const& a)
    : impl_{T(a)}
  {
  }

  template<class T2>
  explicit polynomial(polynomial<T2> const& rhs)
  {
    for (auto const& c : rhs.impl_) {
      impl_.push_back(T(c));
    }
  }
  template<class T2>
  polynomial<T>& operator=(T2 const& rhs)
  {
    impl_ = {T(rhs)};
    return *this;
  }

  template<class T2>
  polynomial<T>& operator=(polynomial<T2> const& rhs)
  {
    impl_.clear();
    for (auto const& c : rhs.impl_) {
      impl_.push_back(T(c));
    }
    return *this;
  }

  std::size_t order() const
  {
    return impl_.size() - 1;
  }

  void order(std::size_t n)
  {
    impl_.resize(n + 1);
  }

  // iterators
  decltype(impl_.begin()) begin() { return impl_.begin(); }
  decltype(impl_.end())   end()   { return impl_.end(); }

  decltype(impl_.cbegin()) begin() const { return impl_.begin(); }
  decltype(impl_.cend())   end()   const { return impl_.end(); }

  // indexing
  T& operator[](std::size_t i)
  {
    return impl_[i];
  }

  const T& operator[](std::size_t i) const
  {
    return impl_[i];
  }

  // comparison
  bool operator==(const polynomial<T>& rhs) const
  {
    return impl_ == rhs.impl_;
  }

  /// arithmetic
  polynomial<T>& operator+=(const T& rhs)
  {
    impl_[0] += rhs;
    return *this;
  }

  polynomial<T>& operator-=(const T& rhs)
  {
    impl_[0] -= rhs;
    return *this;
  }

  polynomial<T>& operator*=(const T& rhs)
  {
    for (auto& c : impl_)
      c *= rhs;
    normalize();
    return *this;
  }

  polynomial<T>& operator+=(const polynomial<T>& rhs)
  {
    order(std::max(order(), rhs.order()));
    for (std::size_t i = 0; i <= rhs.order(); i++)
      impl_[i] += rhs.impl_[i];
    normalize();
    return *this;
  }

  polynomial<T>& operator-=(const polynomial<T>& rhs)
  {
    order(std::max(order(), rhs.order()));
    for (std::size_t i = 0; i <= rhs.order(); i++)
      impl_[i] -= rhs.impl_[i];
    normalize();
    return *this;
  }

  polynomial<T>& operator*=(const polynomial<T>& rhs)
  {
    polynomial<T> product;
    product.order(order() + rhs.order());
    for (std::size_t i = 0; i <= order(); i++) {
      for (std::size_t j = 0; j <= rhs.order(); j++)
        product.impl_[i + j] += impl_[i] * rhs.impl_[j];
    }
    impl_.swap(product.impl_);
    return *this;
  }

  // left_shiftable
  polynomial<T>& operator<<=(std::size_t rhs)
  {
    impl_.insert(impl_.begin(), rhs, T(0));
    return *this;
  }

  const polynomial<T> operator-() const
  {
    polynomial<T> y;
    y.order(order());
    for (std::size_t i = 0; i <= order(); i++)
      y.impl_[i] = -impl_[i];
    return y;
  }

};

template<class T>
std::ostream& operator<<(std::ostream& o, const polynomial<T>& p)
{
  for (auto i = 0u; i <= p.order(); i++) {
    if (p[i] == T(0))
      continue;
    auto c = p[i];
    if (c < 0) {
      o << "- ";
      c = -c;
    } else if (i > 0) {
      o << "+ ";
    }
    if (c != 1 or i == 0) {
      o << c << " ";
    }
    if (i > 0) {
      o << "x";
      if (i > 1)
        o << "^" << i;
      o << " ";
    }
  }
  return o;
}

#endif // POLYNOMIAL_HPP_
