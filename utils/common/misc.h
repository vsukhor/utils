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

/**
* \file misc.h
* A loose collection of functions of common use.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_COMMON_MISC_H
#define UTILS_COMMON_MISC_H

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sys/stat.h>
#include <vector>

//#include <boost/filesystem.hpp>
#include "constants.h"

#ifdef USE_UTILS_XASSERT
    #define XASSERT(EX, msg) \
        (void)((EX) || (Utils::Common::assert_fun(#EX, __FILE__, __LINE__, msg), 0))
#else
    #define XASSERT(EX, msg)
#endif

/// General stuff.
namespace Utils::Common {

/// Convert a double \b number to string.
std::string operator"" _str (long double number);

/// Convert an unsigned long long \b number to string.
std::string operator"" _str (unsigned long long number);

/// Assertion function called from XASSERT macro.
long long assert_fun(
    const char* EX,
    const char *file,
    int line,
    const std::string& msg );

/// Trin the string \b str from whitespaces.
std::string trim(
    const std::string& str,
    const std::string& whitespace  // = " "
);

/// Check that the file named \b name exists.
bool file_exists( const std::string& name );
//bool fileExists( const std::string& name );

/// Check that the directory named \b pathstrconst exists.
bool directory_exists( const std::string& pathstrconst );
//void check_directory(const std::string& s );

/// Cop file named \b fname1 to \b fname2.
void copy_text_file(const std::string& fname1,
                    const std::string& fname2);

/// Sum at compile time.
template <typename T,
          typename Q,
          T (Q::* P)() const> // member function pointer parameter
struct Adder {
    int operator()(const T& i, const Q& o) const {
        return (o.*P)() + i;
    }
};

/// 2D stl-based vectors.
namespace Vec2 {

/// Make a 2D vector having the same dimensions as \b as.
template <typename T1, typename T2> inline
vec2<T1> array_like( const vec2<T2>& as )
{
    vec2<T1> me(as.size());
    for (szt i=0; i<as.size(); i++)
        me[i].resize(as[i].size());
    return me;
}

/// Make a 2D vector having dimensions \b x and \b y eventually initializing it to \b ini.
template <typename T> inline
vec2<T> make( const szt x,
              const szt y,
              const T ini )  // =zero<T>
{
    vec2<T> v(x);
    for (auto& vv : v)
        vv.resize(y, ini);
    return v;
}

/// Return the size of 2D vector \b v.
template <typename T> inline
szt size( const vec2<T>& v )
{
    szt s {};
    for (auto& vv : v)
        s += vv.size();
    return s;
}

/// Add a scalar to 2D vector \b p.
template <typename T>
void add_scalar( const T d,
                 vec2<T>& p )
{
    for (auto& o : p)
        for (auto& oo : o)
            oo += d;
}

/// Fill 2D vector \b v with value \b val.
template <typename T> inline
void fill( vec2<T>& v,
           const T val )
{
    for (auto& o : v)
        for (auto& oo : o)
            oo = val;
}

}  // namespace Vec2

/// 3D stl-based vectors.
namespace Vec3 {

template <typename T> inline
vec3<T> make( const szt x,
              const szt y,
              const szt z,
              const T ini )  // =zero<T>
{
    vec3<T> v(x);
    for (auto& vv : v) {
        vv.resize(y);
        for (auto& vvv : vv)
            vvv.resize(z, ini);
    }
    return v;
}

}  // namespace Vec3

/// \brief Partial (cumulative) sum (in place) of the array \b v.
/// \param u Initial value (of the first element)
/// \param v The array data.
/// \param from First array index to use.
/// \param num Number of array elements to sum up.
template <typename T> inline
void partial_sum( T const* u,
                  T* v,
                  const szt from,
                  const szt num ) noexcept
{
    const auto f {static_cast<int>(from)};
    const auto n {static_cast<int>(num)};

    v[f] = u[f];
    if (f) v[f] += v[f-1];

    for (int i=f+1; i<f+n; i++)
        v[i] = v[i-1] + u[i];
}

/// Average of the vector elements.
template <typename T> inline
auto avg( const std::vector<T>& v )
{
    if constexpr (std::is_integral<T>::value)
        return std::accumulate(v.begin(), v.end(), 0LL) /
               static_cast<double>(v.size());
    else
        return std::accumulate(v.begin(), v.end(), zero<T>) /
               v.size();
}

/// Variance of the vector elements.
template <typename T> inline
T var( const std::vector<T>& n )
{
    T v {};
    T mn {avg(n)};
    for (const auto& o : n)
        v += (o - mn) * (o - mn);

    return v / n.size();
}

/// Average of the vector elements.
template <typename T> constexpr
std::vector<T> exp_num( const T b,
                        const T r,
                        const T dx ) noexcept
{
    uint n {r / dx};
    std::vector<T> x(n);
    std::vector<T> q(n);

    x[0] = zero<T>;
    for (uint i=1; i<n; i++)
        x[i] = x[i-1] + dx;

    const auto c {b / (std::exp(b * r) - one<T>)};

    for (uint i=0; i<n; i++)
        q[i] = c * std::exp(b * x[i]);

    return q;
}

/// \brief Sigmoidal function.
/// \param x Independent variable.
/// \param x0 Position of the turning point.
/// \param steepness Steepness parameter.
template <typename T> constexpr
T sigmoid_decay( const T x,
                 const T x0,
                 const T steepness ) noexcept
{
    return one<T> / (one<T> + std::exp((x - x0) * steepness));
}

/// \brief Find non-zero vector elements.
/// \return how many elements in b /= 0 putting their indices to j
template <typename T> 
szt find( const std::vector<T>& b, 
          std::vector<szt>& j ) noexcept
{
    j.clear();
    for (szt i=0; i<b.size(); i++) 
        if (b[i] != zero<T>)
            j.push_back(i);

    return j.size();
}

/// \brief Find index of the last minimal vector element.
/// \return Index of the last minimal element of vector \b v
template <typename T> constexpr
szt index_min( const std::vector<T>& v ) noexcept
{
    szt k {};
    auto m {v[0]};
    for (szt i=1; i<v.size(); i++)
        if (v[i] < m) {
            k = i;
            m = v[i];
        }
    return k;
}

/// The standard Gaussian function calculated at position \b x.
template <typename T> constexpr
T gaussian( const T x ) noexcept
{
    return std::exp(-x * x / two<T>) /
           std::sqrt(twopi<T>);
}

/// The zero-mean Gaussian function calculated at position \b x.
template <typename T> constexpr
T gaussian( const T x,
            const T var ) noexcept
{
    return std::exp(-x*x / (two<T> * var)) /
           std::sqrt(twopi<T> * var);
}

/// The Gaussian function calculated at position \b x.
template <typename T> constexpr
T gaussian( const T x,
            const T mean,
            const T var ) noexcept
{
    return std::exp(-(x - mean) * (x - mean) / (two<T> * var)) /
           std::sqrt(twopi<T>*var);
}

/// The zero-mean Gaussian function calculated at position \b x.
template <typename T>
T gaussian_fun( T x,
                T mean,
                T sigma )
{
    return one<T> /
           std::sqrt(twopi<T>) *
           std::exp(-(x - mean) * (x - mean) / (sigma * sigma) / two<T>);
}


/// \brief Removes from a vector all instances of an element a.
/// \tparam T Data type (must be EqualityComparable).
/// \param v The vector.
/// \param a The value to remove.
template <typename T>
void remove_vector_element( std::vector<T>& v,
                            const T& a )
{
    for (auto i=v.begin(); i!=v.end(); i++)
        if (a == *i) {
            v.erase(i);
            return;
        }
    XASSERT(false, "Error in remove_vector_element: element not found");
}

}  // namespace Utils::Common


#endif // UTILS_COMMON_MISC_H
