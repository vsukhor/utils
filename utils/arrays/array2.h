/* =============================================================================

 Copyright (C) 2009-2021 Valerii Sukhorukov. All Rights Reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

================================================================================
*/

/// \file array2.h
/// \brief Two-element arrays.
/// \author Valerii Sukhorukov


#ifndef UTILS_ARRAYS_ARRAY2_H
#define UTILS_ARRAYS_ARRAY2_H

#include <array>
#include <cmath>
#include <fstream>
#include <type_traits>
#include <vector>

#include "../common/misc.h"
#include "../constants.h"
#include "_misc.h"

namespace utils::arrays {

/// \brief Two-element arrays.
/// \details This class specializes array template for two-element array of
/// arithmetic types. Implements convenient arithmetics as well as some
/// functionaity commonly used in 2-dimensional geometric applications.
/// \tparam T Type of the elements.
template <typename T>
class array<2,T,std::enable_if_t<std::is_arithmetic_v<T>>> {

static constexpr int len {2};

T n[len] = {};
 
public:

array() noexcept = default;

array( T n ) noexcept
    : n {n, n}
{}

array( T n1, T n2 ) noexcept
    : n {n1, n2}
{}

array( const array& p ) noexcept
    : n {p[0], p[1]}
{}

array( const std::array<T,2>& p ) noexcept
    : n {p[0], p[1]}
{}

template <typename K>
array( const std::array<K,2>& p ) noexcept
    : n {p[0], p[1]}
{}

array( array&& p ) noexcept = default;
array& operator=( array&& p ) noexcept = default;
~array() = default;

template <typename Q>
constexpr array<2,Q> cast_static() const noexcept
{
    return {static_cast<Q>(n[0]),
            static_cast<Q>(n[1])};
}

constexpr array& operator=( const array& p ) noexcept
{
    if (this != &p) {
        n[0] = p[0];
        n[1] = p[1];
    }
    return *this;
}

constexpr array& operator=( const std::array<T,2>& p ) noexcept
{
    if (*this != p) {
        n[0] = p[0];
        n[1] = p[1];
    }
    return *this;
}

constexpr array& operator=( const T p[] ) noexcept
{
    if (n != p) {
        n[0] = p[0];
        n[1] = p[1];
    }
    return *this;
}

constexpr array& operator=( T p ) noexcept
{
    n[0] = p;
    n[1] = p;
    return *this;
}

constexpr array operator+( const array& p ) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1] };
}

template <typename K>
constexpr array operator+( const array<2,K>& p ) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1] };
}

constexpr array operator+( const T p[] ) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1] };
}

constexpr array operator+( T p ) const noexcept
{
    return { n[0] + p,
             n[1] + p };
}

constexpr array& operator+=( const array& p ) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    return *this;
}

constexpr array& operator+=( const T p[] ) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    return *this;
}

constexpr array operator-() const noexcept
{
    return {-n[0], -n[1]};
}

constexpr array operator-( const array& p ) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1] };
}

constexpr array operator-( const T p[] ) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1] };
}

constexpr array& operator-=( const array& p ) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    return *this;
}

constexpr array& operator-=( const T p[] ) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    return *this;
}

constexpr array operator-( T p ) const noexcept
{
    return { n[0] - p,
             n[1] - p };
}

constexpr array operator*( const array& p ) const noexcept
{
    return { n[0] * p[0],
             n[1] * p[1] };
}

constexpr array operator*( const T p[] ) const noexcept
{
    return { n[0]*p[0],
             n[1]*p[1] };
}

constexpr array& operator*=( const array& p ) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    return *this;
}

constexpr array& operator*=( const T p[] ) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    return *this;
}

constexpr array operator*( T p ) const noexcept
{
    return { n[0] * p,
             n[1] * p };
}

constexpr array operator/( const array& p ) const noexcept
{
    return { n[0] / p[0],
             n[1] / p[1] };
}

constexpr array& operator/=( const array& p ) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    return *this;
}

constexpr array& operator/=( const T p[] ) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    return *this;
}

constexpr array operator/( T p ) const noexcept
{
    return { n[0] / p,
             n[1] / p };
}

constexpr bool operator==( const array& p ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1];
}

constexpr bool operator==( const std::array<T,2>& p ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1];
}

constexpr bool operator==( const T p[] ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1];
}

constexpr bool operator==( T p ) const noexcept
{
    return n[0] == p &&
           n[1] == p;
}

constexpr bool operator!=( const array& p ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1];
}

constexpr bool operator!=( const std::array<T,2>& p ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1];
}

constexpr bool operator!=( const T p[] ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1];
}

constexpr bool operator!=( T p ) const noexcept
{
    return n[0] != p ||
           n[1] != p;
}

constexpr bool operator<( const array& p ) const noexcept
{
    return n[0] < p[0] &&
           n[1] < p[1];
}

constexpr bool operator<( T p ) const noexcept
{
    return n[0] < p &&
           n[1] < p;
}

constexpr bool operator<=( const array& p ) const noexcept
{
    return n[0] <= p[0] &&
           n[1] <= p[1];
}

constexpr bool operator<=( T p ) const noexcept
{
    return n[0] <= p &&
           n[1] <= p;
}

constexpr bool operator>( const array& p ) const noexcept
{
    return n[0] > p[0] &&
           n[1] > p[1];
}

constexpr bool operator>( T p ) const noexcept
{
    return n[0] > p &&
           n[1] > p;
}

constexpr bool operator>=( const array& p ) const noexcept
{
    return n[0] >= p[0] &&
           n[1] >= p[1];
}

constexpr bool operator>=( T p ) const noexcept
{
    return n[0] >= p &&
           n[1] >= p;
}

constexpr T operator[]( const int i ) const noexcept
{
    XASSERT(i >= 0 && i < len, "Index out of bounds.");

    return n[i];
}

constexpr T& operator[]( const int i ) noexcept
{
    XASSERT(i >= 0 && i < len, "Index out of bounds.");

    return n[i];
}

[[nodiscard]] constexpr bool contains( T p ) const noexcept
{
    return n[0] == p ||
           n[1] == p;
}

constexpr void reflect() noexcept
{
    T temp = n[0];
    n[0] = n[1];
    n[1] = temp;
}

constexpr int find( T p ) noexcept
{
    return p == n[0] ? 0
                     : (p == n[1] ? 1 : -1);
}

constexpr T other_than( T p ) noexcept
{
    return p == n[0] ? n[1]
                     : (p == n[1] ? n[0] : -1);
}

[[nodiscard]] constexpr T sum() const noexcept
{
    return n[0] + n[1];
}

[[nodiscard]] constexpr T dotpr() const noexcept
{
    return n[0] * n[0] +
           n[1] * n[1];
}

[[nodiscard]] constexpr T dotpr(
    const array& a
) const noexcept
{
    return n[0] * a.n[0] +
           n[1] * a.n[1];
}

static constexpr T dotpr(
    const array& a1,
    const array& a2
) noexcept
{
    return a1.n[0] * a2.n[0] +
           a1.n[1] * a2.n[1];
}

[[nodiscard]] constexpr T norm() const noexcept
{
    return std::sqrt(dotpr());
}

[[nodiscard]] constexpr array unitv() const noexcept
{
    return *this / norm();
}
// Scalar projection of *this onto array b.
[[nodiscard]]constexpr T scaProjection(
    const array& b
) const noexcept
{
    return dotpr(b) / b.norm();
}

// Vector projection of *this onto array b.
[[nodiscard]] constexpr array vecProjection(
    const array& b
) const noexcept
{
    return b.unitv() * scaProjection(b);
}

// Magnitude of the vector that would result from a regular 3D cross product
// of the input vectors, taking their Z values implicitly as 0.
static constexpr T crosspr(
    const array& p1,
    const array& p2
) noexcept
{
    return p1[0] * p2[1] - p1[1] * p2[0];
}

[[nodiscard]] constexpr T crosspr( const array& p ) const noexcept
{
    return n[0] * p[1] - n[1] * p[0];
}

static constexpr bool less0(
    const array& a1,
    const array& a2
) noexcept
{
    return a1[0] < a2[0];
}

static constexpr bool less1(
    const array& a1,
    const array& a2
) noexcept
{
    return a1[1] < a2[1];
}

inline void print(
    std::ostream& os, bool end
) const noexcept
{
    os << n[0] << " " << n[1];
    if (end) os << std::endl;
}

static constexpr array sum(
    const std::vector<array>& a
) noexcept
{
    array res{};
    for (const auto o : a)
        res += o;
    return res;
}

static constexpr array sum(
    const array a[],
    size_t i1,
    size_t i2
) noexcept
{
    array res{};
    for (size_t i=i1; i<=i2; i++)
        res += a[i];
    return res;
}

static constexpr array mean(
    const std::vector<array>& a
) noexcept
{
    return sum(a) / a.size();
}

static constexpr array mean(
    const array a[],
    size_t i1,
    size_t i2
) noexcept
{
    return sum(a, i1, i2) / (i2 - i1 + 1);
}

void read( std::ifstream& ist ) noexcept
{
    ist.read(reinterpret_cast<char*>(&n[0]), sizeof(T));
    ist.read(reinterpret_cast<char*>(&n[1]), sizeof(T));
}
void write( std::ofstream& ost ) const noexcept
{
    ost.write(reinterpret_cast<const char*>(&n[0]), sizeof(T));
    ost.write(reinterpret_cast<const char*>(&n[1]), sizeof(T));
}
};
    
}  // namespace utils::arrays

#endif // UTILS_ARRAYS_ARRAY2_H
