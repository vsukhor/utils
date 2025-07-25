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
 * \file matrix44.h
 * \brief Implementation of 4x4 matrices.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_MATRICES_MATRIX44
#define UTILS_MATRICES_MATRIX44


#include <array>
#include <cmath>
#include <concepts>

#include "../arrays/all.h"
#include "../constants.h"
#include "misc.h"

namespace utils::matrices {


/*
Layout:

00 01 02 03
10 11 12 13
20 21 22 23
30 31 32 33
*/

template<arithmetic ar>
struct Matrix<ar, 4> {

    using A4 = arrays::A4<ar>;

    static constexpr int order {4};

    static constexpr Matrix<ar, order> I
        {{{one<ar>, zero<ar>, zero<ar>, zero<ar>},
          {zero<ar>, one<ar>, zero<ar>, zero<ar>},
          {zero<ar>, zero<ar>, one<ar>, zero<ar>},
          {zero<ar>, zero<ar>, zero<ar>, one<ar>}}};



    constexpr Matrix() = default;
    explicit constexpr Matrix(const ar u[order][order]) noexcept;
    explicit constexpr Matrix(const std::array<A4, order>& a) noexcept;
    explicit constexpr Matrix(const Matrix& m) noexcept;
    explicit constexpr Matrix(Matrix&& m) noexcept;
    constexpr Matrix& operator=(const Matrix& m) const noexcept;
    constexpr Matrix& operator=(Matrix&& m) noexcept;
    ~Matrix() = default;


    constexpr ar& operator()(int i, int j) noexcept;
    constexpr ar operator()(int i, int j) const noexcept;
    constexpr std::array<ar, order> row(int i) const noexcept;
    constexpr void set_row_to(int i,
                              const std::array<ar, order>& r) noexcept;
    constexpr void set_row_to(int i,
                              const ar r[order]) noexcept;
    constexpr std::array<ar, order> col(int j) const noexcept;
    constexpr void set_col_to(int j,
                              const std::array<ar, order>& c) noexcept;
    constexpr void set_col_to(int j,
                              const ar c[order]) noexcept;

    constexpr bool operator==(const Matrix& m) const noexcept;
    constexpr Matrix operator+(const Matrix& m) const noexcept;
    constexpr Matrix operator-(const Matrix& m) const noexcept;
    constexpr Matrix scale(ar c) const noexcept;
    constexpr Matrix ele_mul(const Matrix& m) const noexcept;
    constexpr Matrix mat_mul(const Matrix& m) const noexcept;
    constexpr bool all(std::predicate<ar> auto&& cond) const noexcept;
    constexpr bool any(std::predicate<ar> auto&& cond) const noexcept;

    constexpr ar det() const noexcept;
    constexpr ar det(
        ar a0123,
        ar a0223,
        ar a0323,
        ar a1223,
        ar a1323,
        ar a2323
    ) const noexcept;

    constexpr bool is_singular() const noexcept;
    constexpr bool is_singular(ar det) const noexcept;
    constexpr bool is_orthogonal() const noexcept;
    constexpr bool is_close_to(const Matrix& m,
                               ar tol=EPS<ar>) const noexcept;
    constexpr bool invert(Matrix& inv) const noexcept;
    constexpr Matrix t() const noexcept;
    constexpr std::array<ar, order> diag() const noexcept;
    constexpr ar trace() const noexcept;

private:

    ar u[order][order];   ///< Data. Zero-initialized by default.
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<arithmetic ar>
constexpr
Matrix<ar, 4>::
Matrix(
    const ar _u[order][order]
) noexcept
{
    std::copy(&_u[0][0], &_u[0][0] + order*order, &u[0][0]);
}


template<arithmetic ar>
constexpr
Matrix<ar, 4>::
Matrix(
    const std::array<A4, order>& a
) noexcept
{
    for (int i=0; i<order; i++)
        for (int j=0; j<order; j++)
            u[i][j] = a[i][j];
}


template<arithmetic ar>
constexpr
Matrix<ar, 4>::
Matrix(
    const Matrix& m
) noexcept
{
    if (this != &m)
        std::copy(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);
}


template<arithmetic ar>
constexpr
Matrix<ar, 4>& Matrix<ar, 4>::
operator=(
    const Matrix& m
) const noexcept
{
    if (this != &m)
        std::copy(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);

    return *this;
}


template<arithmetic ar>
constexpr
Matrix<ar, 4>::
Matrix(
    Matrix&& m
) noexcept
{
    std::move(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);
}


template<arithmetic ar>
constexpr
Matrix<ar, 4>& Matrix<ar, 4>::
operator=(
    Matrix&& m
) noexcept
{
    if (this != &m)
        std::move(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);

    return *this;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
operator==(const Matrix& m) const noexcept
{
    return u[0][0] == m.u[0][0] &&
           u[0][1] == m.u[0][1] &&
           u[0][2] == m.u[0][2] &&
           u[0][3] == m.u[0][3] &&
           u[1][0] == m.u[1][0] &&
           u[1][1] == m.u[1][1] &&
           u[1][2] == m.u[1][2] &&
           u[1][3] == m.u[1][3] &&
           u[2][0] == m.u[2][0] &&
           u[2][1] == m.u[2][1] &&
           u[2][2] == m.u[2][2] &&
           u[2][3] == m.u[2][3] &&
           u[3][0] == m.u[3][0] &&
           u[3][1] == m.u[3][1] &&
           u[3][2] == m.u[3][2] &&
           u[3][3] == m.u[3][3];
}


template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
operator+(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] + m.u[0][0],
          u[0][1] + m.u[0][1],
          u[0][2] + m.u[0][2],
          u[0][3] + m.u[0][3]},
         {u[1][0] + m.u[1][0],
          u[1][1] + m.u[1][1],
          u[1][2] + m.u[1][2],
          u[1][3] + m.u[1][3]},
         {u[2][0] + m.u[2][0],
          u[2][1] + m.u[2][1],
          u[2][2] + m.u[2][2],
          u[2][3] + m.u[2][3]},
         {u[3][0] + m.u[3][0],
          u[3][1] + m.u[3][1],
          u[3][2] + m.u[3][2],
          u[3][3] + m.u[3][3]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
operator-(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] - m.u[0][0],
          u[0][1] - m.u[0][1],
          u[0][2] - m.u[0][2],
          u[0][3] - m.u[0][3]},
         {u[1][0] - m.u[1][0],
          u[1][1] - m.u[1][1],
          u[1][2] - m.u[1][2],
          u[1][3] - m.u[1][3]},
         {u[2][0] - m.u[2][0],
          u[2][1] - m.u[2][1],
          u[2][2] - m.u[2][2],
          u[2][3] - m.u[2][3]},
         {u[3][0] - m.u[3][0],
          u[3][1] - m.u[3][1],
          u[3][2] - m.u[3][2],
          u[3][3] - m.u[3][3]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
scale(const ar c) const noexcept
{
    ar k[order][order];
    std::copy(&u[0][0], &u[0][0] + order*order, &k[0][0]);
    for (int i=0; i<order; i++)
        for (int j=0; j<order; j++)
            k[i][j] *= c;

    return Matrix {k};
}

template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
ele_mul(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] * m.u[0][0],
          u[0][1] * m.u[0][1],
          u[0][2] * m.u[0][2],
          u[0][3] * m.u[0][3]},
         {u[1][0] * m.u[1][0],
          u[1][1] * m.u[1][1],
          u[1][2] * m.u[1][2],
          u[1][3] * m.u[1][3]},
         {u[2][0] * m.u[2][0],
          u[2][1] * m.u[2][1],
          u[2][2] * m.u[2][2],
          u[2][3] * m.u[2][3]},
         {u[3][0] * m.u[3][0],
          u[3][1] * m.u[3][1],
          u[3][2] * m.u[3][2],
          u[3][3] * m.u[3][3]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
mat_mul(const Matrix& m) const noexcept
{
    ar k[order][order] {
        {u[0][0]*m.u[0][0] + u[0][1]*m.u[1][0] + u[0][2]*m.u[2][0] + u[0][3]*m.u[3][0],
         u[0][0]*m.u[0][1] + u[0][1]*m.u[1][1] + u[0][2]*m.u[2][1] + u[0][3]*m.u[3][1],
         u[0][0]*m.u[0][2] + u[0][1]*m.u[1][2] + u[0][2]*m.u[2][2] + u[0][3]*m.u[3][2],
         u[0][0]*m.u[0][3] + u[0][1]*m.u[1][3] + u[0][2]*m.u[2][3] + u[0][3]*m.u[3][3]},

        {u[1][0]*m.u[0][0] + u[1][1]*m.u[1][0] + u[1][2]*m.u[2][0] + u[1][3]*m.u[3][0],
         u[1][0]*m.u[0][1] + u[1][1]*m.u[1][1] + u[1][2]*m.u[2][1] + u[1][3]*m.u[3][1],
         u[1][0]*m.u[0][2] + u[1][1]*m.u[1][2] + u[1][2]*m.u[2][2] + u[1][3]*m.u[3][2],
         u[1][0]*m.u[0][3] + u[1][1]*m.u[1][3] + u[1][2]*m.u[2][3] + u[1][3]*m.u[3][3]},

        {u[2][0]*m.u[0][0] + u[2][1]*m.u[1][0] + u[2][2]*m.u[2][0] + u[2][3]*m.u[3][0],
         u[2][0]*m.u[0][1] + u[2][1]*m.u[1][1] + u[2][2]*m.u[2][1] + u[2][3]*m.u[3][1],
         u[2][0]*m.u[0][2] + u[2][1]*m.u[1][2] + u[2][2]*m.u[2][2] + u[2][3]*m.u[3][2],
         u[2][0]*m.u[0][3] + u[2][1]*m.u[1][3] + u[2][2]*m.u[2][3] + u[2][3]*m.u[3][3]},

        {u[3][0]*m.u[0][0] + u[3][1]*m.u[1][0] + u[3][2]*m.u[2][0] + u[3][3]*m.u[3][0],
         u[3][0]*m.u[0][1] + u[3][1]*m.u[1][1] + u[3][2]*m.u[2][1] + u[3][3]*m.u[3][1],
         u[3][0]*m.u[0][2] + u[3][1]*m.u[1][2] + u[3][2]*m.u[2][2] + u[3][3]*m.u[3][2],
         u[3][0]*m.u[0][3] + u[3][1]*m.u[1][3] + u[3][2]*m.u[2][3] + u[3][3]*m.u[3][3]}
        };

    return Matrix {k};
}


template<arithmetic ar>
constexpr
ar& Matrix<ar, 4>::
operator()(const int i, const int j) noexcept
{
    return u[i][j];
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 4>::
operator()(const int i, const int j) const noexcept
{
    return u[i][j];
}


template<arithmetic ar>
constexpr
std::array<ar, 4> Matrix<ar, 4>::
row(const int i) const noexcept
{
    return {{u[i][0], u[i][1], u[i][2], u[i][3]}};
}


template<arithmetic ar>
constexpr
void Matrix<ar, 4>::
set_row_to(const int i,
           const std::array<ar, order>& r) noexcept
{
    u[i][0] = r[0];
    u[i][1] = r[1];
    u[i][2] = r[2];
    u[i][3] = r[3];
}

template<arithmetic ar>
constexpr
void Matrix<ar, 4>::
set_row_to(const int i,
           const ar r[order]) noexcept
{
    u[i][0] = r[0];
    u[i][1] = r[1];
    u[i][2] = r[2];
    u[i][3] = r[3];
}


template<arithmetic ar>
constexpr
std::array<ar, 4> Matrix<ar, 4>::
col(const int j) const noexcept
{
    return {{u[0][j], u[1][j], u[2][j], u[3][j]}};
}


template<arithmetic ar>
constexpr
void Matrix<ar, 4>::
set_col_to(const int j,
           const std::array<ar, order>& c) noexcept
{
    u[0][j] = c[0];
    u[1][j] = c[1];
    u[2][j] = c[2];
    u[3][j] = c[3];
}

template<arithmetic ar>
constexpr
void Matrix<ar, 4>::
set_col_to(const int j,
           const ar c[order]) noexcept
{
    u[0][j] = c[0];
    u[1][j] = c[1];
    u[2][j] = c[2];
    u[3][j] = c[3];
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 4>::
det() const noexcept
{
    const auto a0123 = u[2][0] * u[3][1] - u[2][1] * u[3][0];
    const auto a0223 = u[2][0] * u[3][2] - u[2][2] * u[3][0];
    const auto a0323 = u[2][0] * u[3][3] - u[2][3] * u[3][0];
    const auto a1223 = u[2][1] * u[3][2] - u[2][2] * u[3][1];
    const auto a1323 = u[2][1] * u[3][3] - u[2][3] * u[3][1];
    const auto a2323 = u[2][2] * u[3][3] - u[2][3] * u[3][2];

    return det(a0123, a0223, a0323, a1223, a1323, a2323);
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 4>::
det(
    const ar a0123,
    const ar a0223,
    const ar a0323,
    const ar a1223,
    const ar a1323,
    const ar a2323
) const noexcept
{
    return u[0][0] * (u[1][1] * a2323 - u[1][2] * a1323 + u[1][3] * a1223)
         - u[0][1] * (u[1][0] * a2323 - u[1][2] * a0323 + u[1][3] * a0223)
         + u[0][2] * (u[1][0] * a1323 - u[1][1] * a0323 + u[1][3] * a0123)
         - u[0][3] * (u[1][0] * a1223 - u[1][1] * a0223 + u[1][2] * a0123);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
is_singular() const noexcept
{
    return det() <= EPS<ar>;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
is_singular(const ar det) const noexcept
{
    return det <= EPS<ar>;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
is_close_to(const Matrix<ar, order>& m,
            const ar tol) const noexcept
{
    return (*this - m).all([&](const ar a){ return a <= tol; });
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
is_orthogonal() const noexcept
{
    return is_close_to(this->mat_mul(this->t()), I);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
all(std::predicate<ar> auto&& cond) const noexcept
{
    return cond(u[0][0]) && cond(u[0][1]) && cond(u[0][2]) && cond(u[0][3]) &&
           cond(u[1][0]) && cond(u[1][1]) && cond(u[1][2]) && cond(u[1][3]) &&
           cond(u[2][0]) && cond(u[2][1]) && cond(u[2][2]) && cond(u[2][3]) &&
           cond(u[3][0]) && cond(u[3][1]) && cond(u[3][2]) && cond(u[3][3]);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
any(std::predicate<ar> auto&& cond) const noexcept
{
    return cond(u[0][0]) || cond(u[0][1]) || cond(u[0][2]) || cond(u[0][3]) ||
           cond(u[1][0]) || cond(u[1][1]) || cond(u[1][2]) || cond(u[1][3]) ||
           cond(u[2][0]) || cond(u[2][1]) || cond(u[2][2]) || cond(u[2][3]) ||
           cond(u[3][0]) || cond(u[3][1]) || cond(u[3][2]) || cond(u[3][3]);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 4>::
invert(
    Matrix<ar, order>& inv
) const noexcept
{
    const auto a0123 = u[2][0] * u[3][1] - u[2][1] * u[3][0];
    const auto a0223 = u[2][0] * u[3][2] - u[2][2] * u[3][0];
    const auto a0323 = u[2][0] * u[3][3] - u[2][3] * u[3][0];
    const auto a1223 = u[2][1] * u[3][2] - u[2][2] * u[3][1];
    const auto a1323 = u[2][1] * u[3][3] - u[2][3] * u[3][1];
    const auto a2323 = u[2][2] * u[3][3] - u[2][3] * u[3][2];

    const auto de = det(a0123, a0223, a0323, a1223, a1323, a2323);

    if (is_singular(de))
        return false;

    const auto di = one<ar> / de;

    const auto a0112 = u[1][0] * u[2][1] - u[1][1] * u[2][0];
    const auto a0113 = u[1][0] * u[3][1] - u[1][1] * u[3][0];
    const auto a0212 = u[1][0] * u[2][2] - u[1][2] * u[2][0];
    const auto a0213 = u[1][0] * u[3][2] - u[1][2] * u[3][0];
    const auto a0312 = u[1][0] * u[2][3] - u[1][3] * u[2][0];
    const auto a0313 = u[1][0] * u[3][3] - u[1][3] * u[3][0];
    const auto a1212 = u[1][1] * u[2][2] - u[1][2] * u[2][1];
    const auto a1213 = u[1][1] * u[3][2] - u[1][2] * u[3][1];
    const auto a1312 = u[1][1] * u[2][3] - u[1][3] * u[2][1];
    const auto a1313 = u[1][1] * u[3][3] - u[1][3] * u[3][1];
    const auto a2312 = u[1][2] * u[2][3] - u[1][3] * u[2][2];
    const auto a2313 = u[1][2] * u[3][3] - u[1][3] * u[3][2];

    inv(0,0) = di *  (u[1][1] * a2323 - u[1][2] * a1323 + u[1][3] * a1223);
    inv(0,1) = di * -(u[0][1] * a2323 - u[0][2] * a1323 + u[0][3] * a1223);
    inv(0,2) = di *  (u[0][1] * a2313 - u[0][2] * a1313 + u[0][3] * a1213);
    inv(0,3) = di * -(u[0][1] * a2312 - u[0][2] * a1312 + u[0][3] * a1212);
    inv(1,0) = di * -(u[1][0] * a2323 - u[1][2] * a0323 + u[1][3] * a0223);
    inv(1,1) = di *  (u[0][0] * a2323 - u[0][2] * a0323 + u[0][3] * a0223);
    inv(1,2) = di * -(u[0][0] * a2313 - u[0][2] * a0313 + u[0][3] * a0213);
    inv(1,3) = di *  (u[0][0] * a2312 - u[0][2] * a0312 + u[0][3] * a0212);
    inv(2,0) = di *  (u[1][0] * a1323 - u[1][1] * a0323 + u[1][3] * a0123);
    inv(2,1) = di * -(u[0][0] * a1323 - u[0][1] * a0323 + u[0][3] * a0123);
    inv(2,2) = di *  (u[0][0] * a1313 - u[0][1] * a0313 + u[0][3] * a0113);
    inv(2,3) = di * -(u[0][0] * a1312 - u[0][1] * a0312 + u[0][3] * a0112);
    inv(3,0) = di * -(u[1][0] * a1223 - u[1][1] * a0223 + u[1][2] * a0123);
    inv(3,1) = di *  (u[0][0] * a1223 - u[0][1] * a0223 + u[0][2] * a0123);
    inv(3,2) = di * -(u[0][0] * a1213 - u[0][1] * a0213 + u[0][2] * a0113);
    inv(3,3) = di *  (u[0][0] * a1212 - u[0][1] * a0212 + u[0][2] * a0112);

    return true;
}


template<arithmetic ar>
constexpr
Matrix<ar, 4> Matrix<ar, 4>::
t() const noexcept
{
    ar k[order][order]
        {{u[0][0], u[1][0], u[2][0], u[3][0]},
         {u[0][1], u[1][1], u[2][1], u[3][1]},
         {u[0][2], u[1][2], u[2][2], u[3][2]},
         {u[0][3], u[1][3], u[2][3], u[3][3]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
std::array<ar, 4> Matrix<ar, 4>::
diag() const noexcept
{
    return {{u[0][0], u[1][1], u[2][2], u[3][3]}};
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 4>::
trace() const noexcept
{
    return u[0][0] + u[1][1] + u[2][2] + u[3][3];
}


}  // namespace utils::matrices

#endif  // UTILS_MATRICES_MATRIX44
