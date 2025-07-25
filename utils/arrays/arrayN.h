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
 * \file arrayN.h
 * \brief Arbitrary-size arrays.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_ARRAYS_ARRAYN_H
#define UTILS_ARRAYS_ARRAYN_H

#include <array>
#include <cmath>
#include <fstream>
#include <vector>

#include "../common/misc.h"
#include "../constants.h"
#include "_misc.h"

/// Arbitrary-size arrays.
namespace utils::arrays {

/// Max array length with specific class.
inline constexpr unsigned MAX_LENGTH_SPECIALIZED = 4;


/**
 * \brief Arbitrary-size arrays.
 * \details This class specializes array template for N-element array of
 * arithmetic types. Implements convenient arithmetics as well as some
 * functionaity commonly used in N-dimensional geometric applications.
 * \tparam N Array length.
 * \tparam T Type of the elements.
 */
template<unsigned N,
         arithmetic T>
requires (std::greater_equal<>()(N, MAX_LENGTH_SPECIALIZED+1))
struct array<N, T> {

static constexpr int size {N};

/// Compile-time indexes of the elements.
static constexpr auto ii {make_iota_array<int,N>(0)};

constexpr array() noexcept = default;

/// Scalar constructor.
constexpr array(const T m) noexcept
{
    for (const auto i : ii)
        n[i] = m;
}

constexpr array(const array& p) noexcept
{
    for (const auto i : ii)
        n[i] = p[i];
}

constexpr array(const std::array<T,N>& p) noexcept
{
    for (const auto i : ii)
        n[i] = p[static_cast<std::size_t>(i)];
}

array( array&& p) noexcept = default;
array& operator=( array&& p) noexcept = default;
~array() = default;

constexpr array& operator=(const array& p) noexcept
{
    if (this != &p)
        for (const auto i : ii)
            n[i] = p[i];

    return *this;
}

constexpr array& operator=(const std::array<T,N>& p) noexcept
{
    if (*this != p)
        for (const auto i : ii)
            n[i] = p[i];

    return *this;
}

constexpr array& operator=(const T p[]) noexcept
{
    if (n != p)
        for (const auto i : ii)
            n[i] = p[i];

    return *this;
}

constexpr array& operator=(const T p) noexcept
{
    for (const auto i : ii)
        n[i] = p;

    return *this;
}

constexpr array operator+(const array& p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] + p[i];

    return q;
}

constexpr array operator+(const T p[]) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] + p[i];

    return q;
}

constexpr array& operator+=(const array& p) noexcept
{
    for (const auto i : ii)
        n[i] += p[i];

    return *this;
}

constexpr array& operator+=(const T p[]) noexcept
{
    for (const auto i : ii)
        n[i] += p[i];

    return *this;
}

constexpr array operator+(const T p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] + p;

    return q;
}

constexpr array operator-() const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = -n[i];

    return q;
}

constexpr array operator-(const array& p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] - p[i];

    return q;
}

constexpr array operator-(const T p[]) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] - p[i];

    return q;
}

constexpr array& operator-=(const array& p) noexcept
{
    for (const auto i : ii)
        n[i] -= p[i];

    return *this;
}

constexpr array& operator-=(const T p[]) noexcept
{
    for (const auto i : ii)
        n[i] -= p[i];

    return *this;
}

constexpr array operator-(const T p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] - p;

    return q;
}

constexpr array operator*(const array& p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] * p[i];

    return q;
}

constexpr array operator*(const T p[]) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] * p[i];

    return q;
}

constexpr array& operator*=(const array& p) noexcept
{
    for (const auto i : ii)
        n[i] *= p[i];

    return *this;
}

constexpr array& operator*=(const T p[]) noexcept
{
    for (const auto i : ii)
        n[i] *= p[i];

    return *this;
}

constexpr array operator*(const T p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] * p;

    return q;
}

constexpr array operator/(const array& p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] / p[i];

    return q;
}

constexpr array operator/(const T p[]) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] / p[i];

    return q;
}

constexpr array& operator/=(const array& p) noexcept
{
    for (const auto i : ii)
        n[i] /= p[i];

    return *this;
}

constexpr array& operator/=(const T p[]) noexcept
{
    for (const auto i : ii)
        n[i] /= p[i];

    return *this;
}

constexpr array operator/(const T p) const noexcept
{
    array q;
    for (const auto i : ii)
        q[i] = n[i] / p;

    return q;
}

constexpr bool operator==(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p[i])
            return false;

    return true;
}

constexpr bool operator==(const T p[]) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p[i])
            return false;

    return true;
}

constexpr bool operator==(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p)
            return false;

    return true;
}

constexpr bool operator!=(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p[i])
            return true;

    return false;
}

constexpr bool operator!=(const T p[]) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p[i])
            return true;

    return false;
}

constexpr bool operator!=(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] != p)
            return true;

    return false;
}

constexpr bool operator<(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] >= p[i])
            return false;

    return true;
}

constexpr bool operator<(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] >= p)
            return false;

    return true;
}

constexpr bool operator<=(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] > p[i])
            return false;

    return true;
}

constexpr bool operator<=(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] > p)
            return false;

    return true;
}

constexpr bool operator>(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] <= p[i])
            return false;

    return true;
}

constexpr bool operator>(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] <= p)
            return false;

    return true;
}

constexpr bool operator>=(const array& p) const noexcept
{
    for (const auto i : ii)
        if (n[i] < p[i])
            return false;

    return true;
}

constexpr bool operator>=(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] < p)
            return false;

    return true;
}

constexpr T operator[](const int i) const noexcept
{
    ASSERT(i >= 0 && i < size, "Index out of bounds: ", i);

    return n[i];
}

T& operator[](const int i) noexcept
{
    ASSERT(i >= 0 && i < size, "Index out of bounds: ", i);

    return n[i];
}

constexpr bool contains(const T p) const noexcept
{
    for (const auto i : ii)
        if (n[i] == p)
            return true;

    return false;
}

constexpr T dotpr() const noexcept
{
    T u {};
    for (const auto i : ii)
            u += n[i]*n[i];

    return u;
}

constexpr T dotpr(
    const array& a
) const noexcept
{
    T u {};
    for (const auto i : ii)
            u += n[i]*a[i];

    return u;
}

static constexpr T dotpr(
    const array& a1,
    const array& a2
) noexcept
{
    T u {};
    for (const auto i : ii)
            u += a1[i]*a2[i];

    return u;
}

constexpr T norm() const noexcept
{
    return std::sqrt(dotpr());
}

array unitv() const noexcept
{
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

constexpr T sum() const noexcept
{
    T u {};
    for (const auto i : ii)
            u += n[i];

    return u;
}

void print(
    std::ostream& os,
    const bool end
) const noexcept
{
    for (const auto i : ii)
        os << n[i] << " ";
    if (end) os << std::endl;
}

static constexpr T sum(
    const array a[],
    const size_t i1,
    const size_t i2
) noexcept
{
    T res {};
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

static constexpr T mean(
    const array a[],
    const size_t i1,
    const size_t i2) noexcept
{
    return sum(a, i1, i2) / (i2 - i1 + 1);
}

void read( std::ifstream& ist) noexcept
{
    for (const auto i : ii)
        ist.read(reinterpret_cast<char*>(&n[i]), sizeof(T));
}

void write( std::ofstream& ost) const noexcept
{
    for (const auto i : ii)
        ost.write(reinterpret_cast<const char*>(&n[i]), sizeof(T));
}

private:

T n[size] = {};

};  // struct array<N, T>


/// Input operator.
template<unsigned N,
         arithmetic T>
requires (std::greater_equal<>()(N, MAX_LENGTH_SPECIALIZED+1))
std::istream& operator>>(
    std::istream& is,
    array<N, T>& a
)
{
    for (const auto o: a)
        is >> o;

    return is;
}


/// Output operator.
template<unsigned N,
         arithmetic T>
requires (std::greater_equal<>()(N, MAX_LENGTH_SPECIALIZED+1))
std::ostream& operator<<(
    std::ostream& os,
    const array<N, T>& a
)
{
    a.print(os, false);

    return os;
}

}    // namespace utils::arrays

#endif // UTILS_ARRAYS_ARRAYN_H
