/* =============================================================================

 Copyright (C) 2009-2025 Valerii Sukhorukov. All Rights Reserved.

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

/**
 * \file array4.h
 * \brief Four-element arrays.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_ARRAYS_ARRAY4_H
#define UTILS_ARRAYS_ARRAY4_H

#include <array>
#include <cmath>
#include <fstream>
#include <vector>

#include "../common/misc.h"
#include "../constants.h"
#include "_misc.h"

/// 4-element arrays.
namespace utils::arrays {

/** \brief Four-element arrays.
 * \details This class specializes array template for four-element array of
 * arithmetic types. Implements convenient arithmetics as well as some
 * functionaity commonly used in 4-dimensional geometric applications.
 * \tparam T Type of the elements.
 */
template<arithmetic T>
struct array<4, T> {

using value_type = T;

static constexpr int size {4};

array() noexcept = default;

array(const T m) noexcept
    : n {m, m, m, m}
{}

/// Constructor from explicit values.
array(const T n0, const T n1, const T n2, const T n3) noexcept
    : n {n0, n1, n2, n3}
{}

array(const array<2,T>& n1, const array<2,T>& n2) noexcept
    : n {n1[0], n1[1], n2[0], n2[1]}
{}

array(const T n1, const array<3,T>& n2) noexcept
    : n {n1, n2[0], n2[1], n2[2]}
{}

array(const array<3,T>& n1, T n2) noexcept
    : n {n1[0], n1[1], n1[2], n2}
{}

array(const array& p) noexcept
    : n {p[0], p[1], p[2], p[3]}
{}

array(const std::array<T,4>& p) noexcept
    : n {p[0], p[1], p[2], p[3]}
{}

template<typename Q>
constexpr array<4,Q> cast_static() const noexcept
{
    return {static_cast<Q>(n[0]),
            static_cast<Q>(n[1]),
            static_cast<Q>(n[2]),
            static_cast<Q>(n[3])};
}

array( array&& p) noexcept = default;
array& operator=( array&& p) noexcept = default;
~array() = default;

array& operator=(const array& p) noexcept
{
    if (this != &p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
        n[3] = p[3];
    }

    return *this;
}

constexpr array& operator=(const std::array<T,4>& p) noexcept
{
    if (*this != p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
        n[3] = p[3];
    }
    return *this;
}

constexpr array& operator=(const T p[]) noexcept
{
    if (n != p) {
        n[0] = p[0];
        n[1] = p[1];
        n[2] = p[2];
        n[3] = p[3];
    }
    return *this;
}

constexpr array& operator=(const T p) noexcept
{
    n[0] = p;
    n[1] = p;
    n[2] = p;
    n[3] = p;

    return *this;
}

constexpr array<2,T> operator()(
    const int i1,
    const int i2
) const noexcept
{
    ASSERT(i1 >= 0 && i1 < size, "Index 1 out of bounds: ", i1);
    ASSERT(i2 >= 0 && i2 < size, "Index 2 out of bounds: ", i2);

    return {n[i1], n[i2]};
}

constexpr array<3,T> operator()(
    const int i1,
    const int i2,
    const int i3
) const noexcept
{
    ASSERT(i1 >= 0 && i1 < size, "Index 1 out of bounds: ", i1);
    ASSERT(i2 >= 0 && i2 < size, "Index 2 out of bounds: ", i2);
    ASSERT(i3 >= 0 && i3 < size, "Index 3 out of bounds; ", i3);

    return {n[i1], n[i2], n[i3]};
}

constexpr array operator+(const array& p) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1],
             n[2] + p[2],
             n[3] + p[3] };
}

constexpr array operator+(const T p[]) const noexcept
{
    return { n[0] + p[0],
             n[1] + p[1],
             n[2] + p[2],
             n[3] + p[3] };
}

constexpr array& operator+=(const array& p) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    n[2] += p[2];
    n[3] += p[3];
    return *this;
}

constexpr array& operator+=(const T p[]) noexcept
{
    n[0] += p[0];
    n[1] += p[1];
    n[2] += p[2];
    n[3] += p[3];
    return *this;
}

constexpr array operator+(const T p) const noexcept
{
    return { n[0] + p,
             n[1] + p,
             n[2] + p,
             n[3] + p };
}

constexpr array operator-() const noexcept
{
    array q;
    q[0] = -n[0];
    q[1] = -n[1];
    q[2] = -n[2];
    q[3] = -n[3];
    return q;
}

constexpr array operator-(const array& p) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1],
             n[2] - p[2],
             n[3] - p[3] };
}

constexpr array operator-(const T p[]) const noexcept
{
    return { n[0] - p[0],
             n[1] - p[1],
             n[2] - p[2],
             n[3] - p[3] };
}

constexpr array& operator-=(const array& p) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    n[2] -= p[2];
    n[3] -= p[3];
    return *this;
}

constexpr array& operator-=(const T p[]) noexcept
{
    n[0] -= p[0];
    n[1] -= p[1];
    n[2] -= p[2];
    n[3] -= p[3];
    return *this;
}

constexpr array operator-(const T p) const noexcept
{
    return { n[0] - p,
             n[1] - p,
             n[2] - p,
             n[3] - p };
}

constexpr array operator*(const array& p) const noexcept
{
    return { n[0] * p[0],
             n[1] * p[1],
             n[2] * p[2],
             n[3] * p[3] };
}

constexpr array operator*(const T p[]) const noexcept
{
    return { n[0] * p[0],
             n[1] * p[1],
             n[2] * p[2],
             n[3] * p[3] };
}
constexpr array& operator*=(const array& p) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    n[2] *= p[2];
    n[3] *= p[3];
    return *this;
}

constexpr array& operator*=(const T p[]) noexcept
{
    n[0] *= p[0];
    n[1] *= p[1];
    n[2] *= p[2];
    n[3] *= p[3];
    return *this;
}

constexpr array operator*(const T p) const noexcept
{
    return { n[0] * p,
             n[1] * p,
             n[2] * p,
             n[3] * p };
}

constexpr array operator/(const array& p) const noexcept
{
    return { n[0] / p[0],
             n[1] / p[1],
             n[2] / p[2],
             n[3] / p[3] };
}

constexpr array operator/(const T p[]) const noexcept
{
    return { n[0] / p[0],
             n[1] / p[1],
             n[2] / p[2],
             n[3] / p[3] };
}

constexpr array& operator/=(const array& p) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    n[2] /= p[2];
    n[3] /= p[3];
    return *this;
}

constexpr array& operator/=(const T p[]) noexcept
{
    n[0] /= p[0];
    n[1] /= p[1];
    n[2] /= p[2];
    n[3] /= p[3];
    return *this;
}

constexpr array operator/(const T p) const noexcept
{
    return { n[0] / p,
             n[1] / p,
             n[2] / p,
             n[3] / p };
}

constexpr bool operator==(const array& p) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2] &&
           n[3] == p[3];
}

constexpr bool operator==(const std::array<T,4>& p) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2] &&
           n[3] == p[3];
}

constexpr bool operator==(const T p[]) const noexcept
{
    return n[0] == p[0] &&
           n[1] == p[1] &&
           n[2] == p[2] &&
           n[3] == p[3];
}

constexpr bool operator==(const T p) const noexcept
{
    return n[0] == p &&
           n[1] == p &&
           n[2] == p &&
           n[3] == p;
}

constexpr bool operator!=(const array& p) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1] ||
           n[2] != p[2] ||
           n[3] != p[3];
}

constexpr bool operator!=(const T p[]) const noexcept
{
    return n[0] != p[0] ||
           n[1] != p[1] ||
           n[2] != p[2] ||
           n[3] != p[3];
}

constexpr bool operator!=(const T p) const noexcept
{
    return n[0] != p ||
           n[1] != p ||
           n[2] != p ||
           n[3] != p;
}

constexpr bool operator<(const array& p) const noexcept
{
    return n[0] < p[0] &&
           n[1] < p[1] &&
           n[2] < p[2] &&
           n[3] < p[3];
}

constexpr bool operator<(const T p) const noexcept
{
    return n[0] < p &&
           n[1] < p &&
           n[2] < p &&
           n[3] < p;
}

constexpr bool operator<=(const array& p) const noexcept
{
    return n[0] <= p[0] &&
           n[1] <= p[1] &&
           n[2] <= p[2] &&
           n[3] <= p[3];
}
constexpr bool operator<=(const T p) const noexcept
{
    return n[0] <= p &&
           n[1] <= p &&
           n[2] <= p &&
           n[3] <= p;
}

constexpr bool operator>(const array& p) const noexcept
{
    return n[0] > p[0] &&
           n[1] > p[1] &&
           n[2] > p[2] &&
           n[3] > p[3];
}

constexpr bool operator>(const T p) const noexcept
{
    return n[0] > p &&
           n[1] > p &&
           n[2] > p &&
           n[3] > p;
}

constexpr bool operator>=(const array& p) const noexcept
{
    return n[0] >= p[0] &&
           n[1] >= p[1] &&
           n[2] >= p[2] &&
           n[3] >= p[3];
}

constexpr bool operator>=(const T p) const noexcept
{
    return n[0] >= p &&
           n[1] >= p &&
           n[2] >= p &&
           n[3] >= p;
}

constexpr T operator[](const int i) const noexcept
{
    ASSERT(i >= 0 && i < size, "Index out of bounds: ", i);

    return n[i];
}

constexpr T& operator[](const int i) noexcept
{
    ASSERT(i >= 0 && i < size, "Index out of bounds: ", i);

    return n[i];
}

constexpr bool contains(const T p) const noexcept
{
    return n[0] == p ||
           n[1] == p ||
           n[2] == p ||
           n[3] == p;
}

constexpr void reflect() noexcept
{
    T temp = n[0];
    n[0] = n[3];
    n[3] = temp;
    temp = n[1];
    n[1] = n[2];
    n[2] = temp;
}
constexpr T dotpr() const noexcept
{
    return ( n[0]*n[0] +
             n[1]*n[1] +
             n[2]*n[2] +
             n[3]*n[3] );
}

constexpr T dotpr(const array& a) const noexcept
{
    return ( n[0]*a.n[0] +
             n[1]*a.n[1] +
             n[2]*a.n[2] +
             n[3]*a.n[3] );
}

static constexpr T dotpr(
    const array& a1,
    const array& a2) noexcept
{
    return ( a1.n[0]*a2.n[0] +
             a1.n[1]*a2.n[1] +
             a1.n[2]*a2.n[2] +
             a1.n[3]*a2.n[3] );
}

constexpr T norm() const noexcept
{
    return std::sqrt(dotpr());
}

constexpr array unitv() const noexcept {
    return *this / norm();
}

// Scalar projection of *this onto array b.
constexpr T scaProjection(
    const array& b
) const noexcept
{
    return dotpr(b) / b.norm();
}

// Vector projection of *this onto array b.
constexpr array vecProjection(
    const array& b
) const noexcept
{
    return b.unitv() * scaProjection(b);
}

constexpr int find(const T p) noexcept
{
    return p == n[0] ? 0 :
          (p == n[1] ? 1 :
          (p == n[2] ? 2 :
          (p == n[3] ? 3 : -1)));
}

array<2,int> find(
    const array<2,T>& p
) noexcept
{
    return (p[0] == n[0] && p[1] == n[0])
           ? array<2,int>(0, 1)
           : (p == n[1] ? 1 :
             (p == n[2] ? 2 :
             (p == n[3] ? 3 : -1)));
}

array<3,T> other_than(
    T p
) noexcept
{
    return p == n[0] ? array<3,T>(n[1], n[2], n[3]) :
          (p == n[1] ? array<3,T>(n[0], n[2], n[3]) :
          (p == n[2] ? array<3,T>(n[0], n[1], n[3]) :
          (p == n[3] ? array<3,T>(n[0], n[1], n[2]) :
          array<3,T>(-1))));
}

array<2,T> other_than(
    const array<2,T>& p
) noexcept
{
    return (p == array<2,T>(n[0], n[1]) ||
            p == array<2,T>(n[1], n[0])) ? array<2,T>(n[2], n[3]) :
          ((p == array<2,T>(n[0], n[2]) ||
            p == array<2,T>(n[2], n[0])) ? array<2,T>(n[1], n[3]) :
          ((p == array<2,T>(n[0], n[3]) ||
            p == array<2,T>(n[3], n[0])) ? array<2,T>(n[1], n[2]) :
          ((p == array<2,T>(n[1], n[2]) ||
            p == array<2,T>(n[2], n[1])) ? array<2,T>(n[0], n[3]) :
          ((p == array<2,T>(n[1], n[3]) ||
            p == array<2,T>(n[3], n[1])) ? array<2,T>(n[0], n[2]) :
          ((p == array<2,T>(n[2], n[3]) ||
            p == array<2,T>(n[3], n[2])) ? array<2,T>(n[0], n[1]) :
          array<2,T>(-1))))));
}

T other_than(
    const array<3,T>& p
) noexcept
{
    return (p == array<3,T>(n[0], n[1], n[2]) ||
            p == array<3,T>(n[0], n[2], n[1]) ||
            p == array<3,T>(n[1], n[0], n[2]) ||
            p == array<3,T>(n[1], n[2], n[0]) ||
            p == array<3,T>(n[2], n[0], n[1]) ||
            p == array<3,T>(n[2], n[1], n[0]) ) ? n[3] :
          ((p == array<3,T>(n[0], n[1], n[3]) ||
            p == array<3,T>(n[0], n[3], n[1]) ||
            p == array<3,T>(n[1], n[0], n[3]) ||
            p == array<3,T>(n[1], n[3], n[0]) ||
            p == array<3,T>(n[3], n[0], n[1]) ||
            p == array<3,T>(n[3], n[1], n[0])) ? n[2] :
          ((p == array<3,T>(n[0], n[3], n[2]) ||
            p == array<3,T>(n[0], n[2], n[3]) ||
            p == array<3,T>(n[3], n[0], n[2]) ||
            p == array<3,T>(n[3], n[2], n[0]) ||
            p == array<3,T>(n[2], n[0], n[3]) ||
            p == array<3,T>(n[2], n[3], n[0])) ? n[1] :
          ((p == array<3,T>(n[3], n[1], n[2]) ||
            p == array<3,T>(n[3], n[2], n[1]) ||
            p == array<3,T>(n[1], n[3], n[2]) ||
            p == array<3,T>(n[1], n[2], n[3]) ||
            p == array<3,T>(n[2], n[3], n[1]) ||
            p == array<3,T>(n[2], n[1], n[3])) ? n[0]
            : -constants::one<T>)));
}

constexpr T sum() const noexcept
{
    return n[0] + n[1] + n[2] + n[3];
}

void print(
    std::ostream& os,
    bool end
) const noexcept
{
    os << n[0] << " " << n[1] << " " << n[2] << " " << n[3];
    if (end) os << std::endl;
}

static constexpr array sum(
    const array a[],
    const size_t& i1,
    const size_t& i2
) noexcept
{
    array res{};
    for (size_t i=i1; i<=i2; i++)
        res += a[i];

    return res;
}

constexpr static array mean(
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

void read(
    std::ifstream& ist
) noexcept
{
    ist.read( reinterpret_cast<char*>(&n[0]), sizeof(T));
    ist.read( reinterpret_cast<char*>(&n[1]), sizeof(T));
    ist.read( reinterpret_cast<char*>(&n[2]), sizeof(T));
    ist.read( reinterpret_cast<char*>(&n[3]), sizeof(T));
}

void write(
    std::ofstream& ost
) const noexcept
{
    ost.write( reinterpret_cast<const char*>(&n[0]), sizeof(T));
    ost.write( reinterpret_cast<const char*>(&n[1]), sizeof(T));
    ost.write( reinterpret_cast<const char*>(&n[2]), sizeof(T));
    ost.write( reinterpret_cast<const char*>(&n[3]), sizeof(T));
}

private:

T n[size] = {};

};  // struct array<4, T>


/// Input operator.
template<arithmetic T>
std::istream& operator>>(
    std::istream& is,
    array<4, T>& a
)
{
    is >> a[0] >> a[1] >> a[2] >> a[3];

    return is;
}


/// Output operator.
template<arithmetic T>
std::ostream& operator<<(
    std::ostream& os,
    const array<4, T>& a
)
{
    a.print(os, false);

    return os;
}


}  // namespace utils::arrays

#endif  // UTILS_ARRAYS_ARRAY4_H
