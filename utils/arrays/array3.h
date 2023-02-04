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

/// \file array3.h
/// \brief Three-element arrays.
/// \author Valerii Sukhorukov


#ifndef UTILS_ARRAYS_ARRAY3_H
#define UTILS_ARRAYS_ARRAY3_H

#include <array>
#include <cmath>
#include <fstream>
#include <type_traits>
#include <vector>

#include "../common/misc.h"
#include "../constants.h"
#include "_misc.h"
#include "array2.h"

/// Library-wide.
namespace utils {

namespace common {

template<typename>
struct Geometric;

}  // namespace common

/// 3-element arrays.
namespace arrays {

/// \brief Three-element arrays.
/// \details This class specializes array template for three-element array
/// of arithmetic types. Implements convenient arithmetics as well as some
/// functionaity commonly used in 3-dimensional geometric applications.
/// \tparam T Type of the elements.
template<typename T>
class array<3,T,std::enable_if_t<std::is_arithmetic<T>::value>> {

static constexpr int len {3};

T n[len] = {};
 
public:

array() noexcept = default;

constexpr
array( const T m ) noexcept
    : n {m, m, m}
{}

constexpr
array( const T n1, const T n2, const T n3 ) noexcept
    : n {n1, n2, n3}
{}

constexpr
array( const T n1, const array<2,T>& n2 ) noexcept
    : n {n1, n2[0], n2[1]}
{}

constexpr
array( const array<2,T>& n1, T n2 ) noexcept
    : n {n1[0], n1[1], n2}
{}

constexpr
array( const array& p ) noexcept
    : n {p[0], p[1], p[2]}
{}

constexpr
array( const std::array<T,3>& p ) noexcept
    : n {p[0], p[1], p[2]}
{}

array( array&& p ) noexcept = default;

array& operator=( array&& p ) noexcept = default;

~array() = default;

explicit array( std::ifstream& ist )
{
    read(ist);
}

template<typename Q>
constexpr
array<3,Q> cast_static() const noexcept
{
    return {static_cast<Q>(n[0]),
            static_cast<Q>(n[1]),
            static_cast<Q>(n[2])};
}

constexpr
array<2,T> operator()(
    const int i1,
    const int i2
) const noexcept
{
    XASSERT(i1 >= 0 && i1 < len, "Index 1 out of bounds.");
    XASSERT(i2 >= 0 && i2 < len, "Index 2 out of bounds.");
    
    return {n[i1], n[i2]};
}

constexpr
T const* data() const noexcept
{
    return n;
}

constexpr
array& operator=( const array& p ) noexcept
{
    if (this != &p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
    }
    return *this;
}

constexpr
array& operator=( const std::array<T,3>& p ) noexcept
{
    if (*this != p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
    }
    return *this;
}

constexpr
array& operator=( const T p[] ) noexcept
{
    if (n != p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
    }
    return *this;
}

constexpr
array& operator=( const T p ) noexcept
{
    n[0] = p;
    n[1] = p;
    n[2] = p;
    
    return *this;
}

constexpr
array operator+( const array& p ) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1],
             n[2] + p[2] };
}

constexpr
array operator+( const T p[] ) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1],
             n[2] + p[2] };
}

constexpr
array& operator+=( const array& p ) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    n[2] += p[2];
    return *this;
}

constexpr
array& operator+=( const T p[] ) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    n[2] += p[2];
    return *this;
}

constexpr
array operator+( const T p ) const noexcept
{
    return { n[0] + p,
             n[1] + p,
             n[2] + p };
}

constexpr
array operator-() const noexcept
{
    return {-n[0], -n[1], -n[2]};
}

constexpr
array operator-( const array& p ) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1],
             n[2] - p[2] };
}

constexpr
array operator-( const T p[] ) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1],
             n[2] - p[2] };
}

constexpr
array& operator-=( const array& p ) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    n[2] -= p[2];
    return *this;
}

constexpr array& operator-=( const T p[] ) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    n[2] -= p[2];
    return *this;
}

constexpr
array operator-( const T p ) const noexcept
{
    return { n[0] - p,
             n[1] - p,
             n[2] - p };
}

constexpr
array operator*( const array& p ) const noexcept
{
    return { n[0] * p[0],
             n[1] * p[1],
             n[2] * p[2] };
}

constexpr
array operator*( const T p[] ) const noexcept
{
    return { n[0] * p[0],
             n[1] * p[1],
             n[2] * p[2] };
}

constexpr
array& operator*=( const array& p ) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    n[2] *= p[2];
    return *this;
}

constexpr
array& operator*=( const T p[] ) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    n[2] *= p[2];
    return *this;
}

constexpr
array operator*( const T p ) const noexcept
{
    return { n[0] * p,
             n[1] * p,
             n[2] * p };
}

constexpr
array operator/( const array& p ) const noexcept
{
    return { n[0] / p[0],
             n[1] / p[1],
             n[2] / p[2] };
}

constexpr
array operator/( const T p[] ) const noexcept
{
    return { n[0] / p[0],
             n[1] / p[1],
             n[2] / p[2] };
}

constexpr
array& operator/=( const array& p ) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    n[2] /= p[2];
    return *this;
}

constexpr
array& operator/=( const T p[] ) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    n[2] /= p[2];
    return *this;
}

constexpr
array operator/( const T p ) const noexcept
{
    return { n[0] / p,
             n[1] / p,
             n[2] / p };
}

constexpr
bool operator==( const array& p ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2];
}

constexpr
bool operator==( const std::array<T,3>& p ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2];
}

constexpr
bool operator==( const T p[] ) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2];
}

constexpr
bool operator==( const T p ) const noexcept
{
    return n[0] == p &&
           n[1] == p &&
           n[2] == p;
}

constexpr
bool operator!=( const array& p ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1] ||
           n[2] != p[2];
}

constexpr
bool operator!=( const std::array<T,3>& p ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1] ||
           n[2] != p[2];
}

constexpr
bool operator!=( const T p[] ) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1] ||
           n[2] != p[2];
}

constexpr
bool operator!=( const T p ) const noexcept
{
    return n[0] != p ||
           n[1] != p ||
           n[2] != p;
}

constexpr
bool operator<( const array& p ) const noexcept
{
    return n[0] < p[0] &&
           n[1] < p[1] &&
           n[2] < p[2];
}

constexpr
bool operator<( const T p ) const noexcept
{
    return n[0] < p &&
           n[1] < p &&
           n[2] < p;
}

constexpr
bool operator<=( const array& p ) const noexcept
{
    return n[0] <= p[0] &&
           n[1] <= p[1] &&
           n[2] <= p[2];
}

constexpr
bool operator<=( const T p ) const noexcept
{
    return n[0] <= p &&
           n[1] <= p &&
           n[2] <= p;
}

constexpr
bool operator>( const array& p ) const noexcept
{
    return n[0] > p[0] &&
           n[1] > p[1] &&
           n[2] > p[2];
}

constexpr
bool operator>( const T p ) const noexcept
{
    return n[0] > p &&
           n[1] > p &&
           n[2] > p;
}

constexpr
bool operator>=( const array& p ) const noexcept
{
    return n[0] >= p[0] &&
           n[1] >= p[1] &&
           n[2] >= p[2];
}

constexpr
bool operator>=( const T p ) const noexcept
{
    return n[0] >= p &&
           n[1] >= p &&
           n[2] >= p;
}

constexpr
T operator[]( const int i ) const noexcept
{
    XASSERT(i >= 0 && i < len, "Index out of bounds.");

    return n[i];
}

T& operator[]( const int i ) noexcept
{
    XASSERT(i >= 0 && i < len, "Index out of bounds.");

    return n[i];
}

const T* at( const int i ) const noexcept
{
    XASSERT(i >= 0 && i < len, "Index out of bounds.");

    return &n[i];
}

constexpr
bool contains( const T p ) const noexcept
{
    return n[0] == p ||
           n[1] == p ||
           n[2] == p;
}

constexpr
void reflect() noexcept
{
    T temp = n[0];
    n[0] = n[2];
    n[2] = temp;
}

constexpr
T dotpr() const noexcept
{
    return ( n[0] * n[0] +
             n[1] * n[1] +
             n[2] * n[2] );
}

constexpr
T dotpr( const array& a ) const noexcept
{
    return ( n[0] * a.n[0] +
             n[1] * a.n[1] +
             n[2] * a.n[2] );
}

static constexpr
T dotpr(
        const array& a1,
        const array& a2
) noexcept
{
    return ( a1.n[0] * a2.n[0] +
             a1.n[1] * a2.n[1] +
             a1.n[2] * a2.n[2] );
}

constexpr
T norm() const noexcept
{
    return std::sqrt(dotpr());
}

constexpr
array unitv() const noexcept
{
    return *this / norm();
}

// Scalar projection of *this onto array b.
constexpr
T scaProjection(
    const array& b
) const noexcept
{
    return dotpr(b) / b.norm();
}

// Vector projection of *this onto array b.
constexpr
array vecProjection(
    const array& b
) const noexcept
{
    return b.unitv() * scaProjection(b);
}

constexpr
int find(
    const T p
) noexcept
{
    return p == n[0] ? 0 :
          (p == n[1] ? 1 :
          (p == n[2] ? 2 : -1));
}

array<2,T> other_than(
    const T p
) noexcept
{
    return p == n[0] ? array<2,T>(n[1], n[2]) :
          (p == n[1] ? array<2,T>(n[0], n[2]) :
          (p == n[2] ? array<2,T>(n[0], n[1]) : array<2,T>(-1)));
}

T other_than( const array<2,T>& p
) noexcept
{
    return (p == array<2,T>(n[0], n[1]) ||
            p == array<2,T>(n[1], n[0])) ? n[2] :
          ((p == array<2,T>(n[0], n[2]) ||
            p == array<2,T>(n[2], n[0])) ? n[1] :
          ((p == array<2,T>(n[1], n[2]) ||
            p == array<2,T>(n[2], n[1])) ? n[0] :
            -constants::one<T>));
}

constexpr
T sum() const noexcept
{
    return n[0] + n[1] + n[2];
}

static constexpr
array crosspr(
    const array& p1,
    const array& p2
) noexcept
{
    return { p1[1] * p2[2] - p1[2] * p2[1],
             p1[2] * p2[0] - p1[0] * p2[2],
             p1[0] * p2[1] - p1[1] * p2[0] };
}

constexpr
array crosspr(
    const array& p2
) const noexcept
{
    return { n[1] * p2[2] - n[2] * p2[1],
             n[2] * p2[0] - n[0] * p2[2],
             n[0] * p2[1] - n[1] * p2[0] };
}

constexpr
void rotate(
    const T rm[3][3]
) noexcept
{
//    const array m {rm[0][0]*n[0]+rm[0][1]*n[1]+rm[0][2]*n[2],
//                    rm[1][0]*n[0]+rm[1][1]*n[1]+rm[1][2]*n[2],
//                    rm[2][0]*n[0]+rm[2][1]*n[1]+rm[2][2]*n[2] }; new
    const array m {rm[0][0] * n[0] + rm[1][0] * n[1] + rm[2][0] * n[2],
                   rm[0][1] * n[0] + rm[1][1] * n[1] + rm[2][1] * n[2],
                   rm[0][2] * n[0] + rm[1][2] * n[1] + rm[2][2] * n[2] };
    *this = m;
}

// Rotate vector v with rotation matrix rm.
static constexpr
array rotate(
    const array& v,
    const T rm[3][3]
) noexcept
{
    return { rm[0][0] * v[0] + rm[1][0] * v[1] + rm[2][0] * v[2],
             rm[0][1] * v[0] + rm[1][1] * v[1] + rm[2][1] * v[2],
             rm[0][2] * v[0] + rm[1][2] * v[1] + rm[2][2] * v[2] };
}

// Rotate vector v through an 'angle' around an 'axis' in 3D.
static constexpr
array rotate(
    const array& axis,
    const array& v,
    const T angle
) noexcept
{
    T rm[3][3];
    common::Geometric<T>::rotmat(axis, angle, rm);

    return { rm[0][0] * v[0] + rm[1][0] * v[1] + rm[2][0] * v[2],
             rm[0][1] * v[0] + rm[1][1] * v[1] + rm[2][1] * v[2],
             rm[0][2] * v[0] + rm[1][2] * v[1] + rm[2][2] * v[2] };
}

constexpr array
rotate(
    const T rm[][3]
) const noexcept
{
    return { rm[0][0] * n[0] + rm[1][0] * n[1] + rm[2][0] * n[2],
             rm[0][1] * n[0] + rm[1][1] * n[1] + rm[2][1] * n[2],
             rm[0][2] * n[0] + rm[1][2] * n[1] + rm[2][2] * n[2] };
}

constexpr
array rotate(
    const array& axis,
    const T angle
) const noexcept
{
    T rm[3][3];
    common::Geometric<T>::rotmat(axis, angle, rm);

    return { rm[0][0] * n[0] + rm[1][0] * n[1] + rm[2][0] * n[2],
             rm[0][1] * n[0] + rm[1][1] * n[1] + rm[2][1] * n[2],
             rm[0][2] * n[0] + rm[1][2] * n[1] + rm[2][2] * n[2] };
}

void print(
    std::ostream& os,
    const bool end
) const noexcept
{
    os << n[0] << " " << n[1] << " " << n[2];
    if (end) os << std::endl;
}

static array sum(
    const array a[],
    const size_t i1,
    const size_t i2
) noexcept
{
    array res {};
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
    const size_t& i1,
    const size_t& i2
) noexcept
{
    return sum(a, i1, i2) / (i2 - i1 + 1);
}

void read( std::ifstream& ist ) noexcept
{
    ist.read(reinterpret_cast<char*>(&n[0]), sizeof(T));
    ist.read(reinterpret_cast<char*>(&n[1]), sizeof(T));
    ist.read(reinterpret_cast<char*>(&n[2]), sizeof(T));
}

void write( std::ofstream& ost ) const noexcept
{
    ost.write(reinterpret_cast<const char*>(&n[0]), sizeof(T));
    ost.write(reinterpret_cast<const char*>(&n[1]), sizeof(T));
    ost.write(reinterpret_cast<const char*>(&n[2]), sizeof(T));
}

[[nodiscard]] constexpr
array<2,int> equal_dims() const noexcept
{
    if (n[0] == n[1]) return {0, 1};
    if (n[0] == n[2]) return {0, 2};
    if (n[1] == n[2]) return {1, 2};

    return {-1, -1};
}

[[nodiscard]] constexpr
int index_min() const noexcept
{
    return (n[0] > n[1]) ? ((n[1] > n[2]) ? 2 : 1)
                         : ((n[0] > n[2]) ? 2 : 0);
}


};

}    // namespace arrays
}    // namespace utils

#endif // UTILS_ARRAYS_ARRAY3_H
