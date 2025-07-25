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
 * \file matrix33.h
 * \brief Implementation of 3x3 matrices.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_MATRICES_MATRIX3
#define UTILS_MATRICES_MATRIX3


#include <array>
#include <cmath>
#include <concepts>

#include "../arrays/all.h"
#include "../constants.h"
#include "misc.h"


namespace utils::matrices {


/*
Layout:

00 01 02
10 11 12
20 21 22
*/

template<arithmetic ar>
struct Matrix<ar, 3> {

    using A3 = arrays::A3<ar>;

    static constexpr int order {3};

    static constexpr Matrix<ar, order> I
        {{{one<ar>, zero<ar>, zero<ar>},
          {zero<ar>, one<ar>, zero<ar>},
          {zero<ar>, zero<ar>, one<ar>}}};


    constexpr Matrix() = default;
    explicit constexpr Matrix(const ar u[order][order]) noexcept;
    explicit constexpr Matrix(const std::array<A3, order>& a) noexcept;
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
        ar a00,
        ar a10,
        ar a20
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
Matrix<ar, 3>::
Matrix(
    const ar _u[order][order]
) noexcept
{
    std::copy(&_u[0][0], &_u[0][0] + order*order, &u[0][0]);
}

template<arithmetic ar>
constexpr
Matrix<ar, 3>::
Matrix(
    const Matrix& m
) noexcept
{
    if (this != &m)
        std::copy(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);
}


template<arithmetic ar>
constexpr
Matrix<ar, 3>::
Matrix(
    const std::array<A3, order>& a
) noexcept
{
    for (int i=0; i<order; i++)
        for (int j=0; j<order; j++)
            u[i][j] = a[i][j];
}


template<arithmetic ar>
constexpr
Matrix<ar, 3>& Matrix<ar, 3>::
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
Matrix<ar, 3>::
Matrix(
    Matrix&& m
) noexcept
{
    std::move(&m.u[0][0], &m.u[0][0] + order*order, &u[0][0]);
}


template<arithmetic ar>
constexpr
Matrix<ar, 3>& Matrix<ar, 3>::
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
bool Matrix<ar, 3>::
operator==(const Matrix& m) const noexcept
{
    return u[0][0] == m.u[0][0] &&
           u[0][1] == m.u[0][1] &&
           u[0][2] == m.u[0][2] &&
           u[1][0] == m.u[1][0] &&
           u[1][1] == m.u[1][1] &&
           u[1][2] == m.u[1][2] &&
           u[2][0] == m.u[2][0] &&
           u[2][1] == m.u[2][1] &&
           u[2][2] == m.u[2][2];
}


template<arithmetic ar>
constexpr
Matrix<ar, 3> Matrix<ar, 3>::
operator+(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] + m.u[0][0],
          u[0][1] + m.u[0][1],
          u[0][2] + m.u[0][2]},
         {u[1][0] + m.u[1][0],
          u[1][1] + m.u[1][1],
          u[1][2] + m.u[1][2]},
         {u[2][0] + m.u[2][0],
          u[2][1] + m.u[2][1],
          u[2][2] + m.u[2][2]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 3> Matrix<ar, 3>::
operator-(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] - m.u[0][0],
          u[0][1] - m.u[0][1],
          u[0][2] - m.u[0][2]},
         {u[1][0] - m.u[1][0],
          u[1][1] - m.u[1][1],
          u[1][2] - m.u[1][2]},
         {u[2][0] - m.u[2][0],
          u[2][1] - m.u[2][1],
          u[2][2] - m.u[2][2]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 3> Matrix<ar, 3>::
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
Matrix<ar, 3> Matrix<ar, 3>::
ele_mul(const Matrix& m) const noexcept
{
    ar k[order][order]
        {{u[0][0] * m.u[0][0],
          u[0][1] * m.u[0][1],
          u[0][2] * m.u[0][2]},
         {u[1][0] * m.u[1][0],
          u[1][1] * m.u[1][1],
          u[1][2] * m.u[1][2]},
         {u[2][0] * m.u[2][0],
          u[2][1] * m.u[2][1],
          u[2][2] * m.u[2][2]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
Matrix<ar, 3> Matrix<ar, 3>::
mat_mul(const Matrix& m) const noexcept
{
    ar k[order][order] {
        {u[0][0]*m.u[0][0] + u[0][1]*m.u[1][0] + u[0][2]*m.u[2][0],
         u[0][0]*m.u[0][1] + u[0][1]*m.u[1][1] + u[0][2]*m.u[2][1],
         u[0][0]*m.u[0][2] + u[0][1]*m.u[1][2] + u[0][2]*m.u[2][2]},

        {u[1][0]*m.u[0][0] + u[1][1]*m.u[1][0] + u[1][2]*m.u[2][0],
         u[1][0]*m.u[0][1] + u[1][1]*m.u[1][1] + u[1][2]*m.u[2][1],
         u[1][0]*m.u[0][2] + u[1][1]*m.u[1][2] + u[1][2]*m.u[2][2]},

        {u[2][0]*m.u[0][0] + u[2][1]*m.u[1][0] + u[2][2]*m.u[2][0],
         u[2][0]*m.u[0][1] + u[2][1]*m.u[1][1] + u[2][2]*m.u[2][1],
         u[2][0]*m.u[0][2] + u[2][1]*m.u[1][2] + u[2][2]*m.u[2][2]}
    };

    return Matrix {k};
}


template<arithmetic ar>
constexpr
ar& Matrix<ar, 3>::
operator()(const int i, const int j) noexcept
{
    return u[i][j];
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 3>::
operator()(const int i, const int j) const noexcept
{
    return u[i][j];
}


template<arithmetic ar>
constexpr
std::array<ar, 3> Matrix<ar, 3>::
row(const int i) const noexcept
{
    return {{u[i][0], u[i][1], u[i][2]}};
}


template<arithmetic ar>
constexpr
void Matrix<ar, 3>::
set_row_to(const int i,
           const std::array<ar, order>& r) noexcept
{
    u[i][0] = r[0];
    u[i][1] = r[1];
    u[i][2] = r[2];
}

template<arithmetic ar>
constexpr
void Matrix<ar, 3>::
set_row_to(const int i,
           const ar r[order]) noexcept
{
    u[i][0] = r[0];
    u[i][1] = r[1];
    u[i][2] = r[2];
}


template<arithmetic ar>
constexpr
std::array<ar, 3> Matrix<ar, 3>::
col(const int j) const noexcept
{
    return {{u[0][j], u[1][j], u[2][j]}};
}


template<arithmetic ar>
constexpr
void Matrix<ar, 3>::
set_col_to(const int j,
           const std::array<ar, order>& c) noexcept
{
    u[0][j] = c[0];
    u[1][j] = c[1];
    u[2][j] = c[2];
}


template<arithmetic ar>
constexpr
void Matrix<ar, 3>::
set_col_to(const int j,
           const ar c[order]) noexcept
{
    u[0][j] = c[0];
    u[1][j] = c[1];
    u[2][j] = c[2];
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 3>::
det() const noexcept
{
    const auto a00 =  (u[1][1]*u[2][2] - u[2][1]*u[1][2]);
    const auto a10 = -(u[1][0]*u[2][2] - u[2][0]*u[1][2]);
    const auto a20 =  (u[1][0]*u[2][1] - u[2][0]*u[1][1]);

    return det(a00, a10, a20);
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 3>::
det(
    const ar a00,
    const ar a10,
    const ar a20
) const noexcept
{
    return u[0][0] * a00 + u[0][1] * a10 + u[0][2] * a20;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
is_singular() const noexcept
{
    return det() <= EPS<ar>;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
is_singular(const ar det) const noexcept
{
    return det <= EPS<ar>;
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
is_close_to(const Matrix<ar, order>& m,
            const ar tol) const noexcept
{
    return (*this - m).all([&](const ar a){ return a <= tol; });
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
is_orthogonal() const noexcept
{
    return is_close_to(this->mat_mul(this->t()), I);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
all(std::predicate<ar> auto&& cond) const noexcept
{
    return cond(u[0][0]) && cond(u[0][1]) && cond(u[0][2]) &&
           cond(u[1][0]) && cond(u[1][1]) && cond(u[1][2]) &&
           cond(u[2][0]) && cond(u[2][1]) && cond(u[2][2]);
}


template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
any(std::predicate<ar> auto&& cond) const noexcept
{
    return cond(u[0][0]) || cond(u[0][1]) || cond(u[0][2]) ||
           cond(u[1][0]) || cond(u[1][1]) || cond(u[1][2]) ||
           cond(u[2][0]) || cond(u[2][1]) || cond(u[2][2]);
}

template<arithmetic ar>
constexpr
bool Matrix<ar, 3>::
invert(
    Matrix<ar, order>& inv
) const noexcept
{
    const auto a00 =  (u[1][1]*u[2][2] - u[2][1]*u[1][2]);
    const auto a10 = -(u[1][0]*u[2][2] - u[2][0]*u[1][2]);
    const auto a20 =  (u[1][0]*u[2][1] - u[2][0]*u[1][1]);

    const auto de = det(a00, a10, a20);

    if (is_singular(de))
        return false;

    const auto di = one<ar> / de;

    const auto a01 = -(u[0][1]*u[2][2] - u[2][1]*u[0][2]);
    const auto a11 =  (u[0][0]*u[2][2] - u[2][0]*u[0][2]);
    const auto a21 = -(u[0][0]*u[2][1] - u[2][0]*u[0][1]);
    const auto a02 =  (u[0][1]*u[1][2] - u[1][1]*u[0][2]);
    const auto a12 = -(u[0][0]*u[1][2] - u[1][0]*u[0][2]);
    const auto a22 =  (u[0][0]*u[1][1] - u[1][0]*u[0][1]);

    inv.u[0][0] = di * a00;
    inv.u[0][1] = di * a01;
    inv.u[0][2] = di * a02;
    inv.u[1][0] = di * a10;
    inv.u[1][1] = di * a11;
    inv.u[1][2] = di * a12;
    inv.u[2][0] = di * a20;
    inv.u[2][1] = di * a21;
    inv.u[2][2] = di * a22;

    return true;
}


template<arithmetic ar>
constexpr
Matrix<ar, 3> Matrix<ar, 3>::
t() const noexcept
{
    ar k[order][order]
        {{u[0][0], u[1][0], u[2][0]},
         {u[0][1], u[1][1], u[2][1]},
         {u[0][2], u[1][2], u[2][2]}};

    return Matrix {k};
}


template<arithmetic ar>
constexpr
std::array<ar, 3> Matrix<ar, 3>::
diag() const noexcept
{
    return {{u[0][0], u[1][1], u[2][2]}};
}


template<arithmetic ar>
constexpr
ar Matrix<ar, 3>::
trace() const noexcept
{
    return u[0][0] + u[1][1] + u[2][2];
}


}  // namespace utils::matrices

#endif  // UTILS_MATRICES_MATRIX3
