/* ==============================================================================

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

============================================================================== */

/**
* \file constants.h
* \brief A collection of useful definitions.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_COMMON_CONSTANTS_H
#define UTILS_COMMON_CONSTANTS_H

#include <limits>
#include <vector>
#include <array>
#include <memory>
#include <numeric>
#include <cmath>

/// Library-wide.
namespace Utils {
/// General stuff.
namespace Common {

extern const std::string SLASH;

#define STR(x) std::to_string(x)

using ulong = unsigned long;
using uint = unsigned int;
using szt = std::size_t;

// container typedefs ==========================================================

// std vector-based 2, 3, 4-dim containers
template <typename T> using vec2 = std::vector<std::vector<T>>;
template <typename T> using vec3 = std::vector<vec2<T>>;
template <typename T> using vec4 = std::vector<vec3<T>>;

// std vector and std array
template <typename T, auto N> using arrvec = std::array<std::vector<T>,N>;
template <typename T, auto N> using vecarr = std::vector<std::array<T,N>>;

// std vector of unique pointers
template <typename T> using vup = std::vector<std::unique_ptr<T>>;

// common constants ============================================================

template <typename T> constexpr T zero {static_cast<T>(0.L)};
template <typename T> constexpr T half {static_cast<T>(.5L)};
template <typename T> constexpr T thrd {static_cast<T>(1.L/3.L)};
template <typename T> constexpr T one {static_cast<T>(1.L)};
template <typename T> constexpr T two {static_cast<T>(2.L)};
template <typename T> constexpr T three {static_cast<T>(3.L)};
template <typename T> constexpr T four {static_cast<T>(4.L)};
template <typename T> constexpr T five {static_cast<T>(5.L)};
template <typename T> constexpr T six {static_cast<T>(6.L)};
template <typename T> constexpr T ten {static_cast<T>(10.L)};

template <typename T> constexpr T
    pi {static_cast<T>(3.1415926535897932384626433832795L)};
template <typename T> constexpr T
    twopi {two<T>*pi<T>};
template <typename T> constexpr T
    halfpi {half<T>*pi<T>};
template <typename T> constexpr T
    sqrtPI {static_cast<T>(1.7724538509055160272981674833411L)};

template <typename T, typename Enabler = void> constexpr T sqrt2PI;
template <typename T> constexpr T
    sqrt2PI<T,std::enable_if_t<std::is_arithmetic<T>::value>> {
        std::pow(twopi<T>, half<T>)
    };

// Templates for numerical limits. =============================================

template <typename T, typename Enabler = void> constexpr T EPS;
template <typename T> constexpr T
    EPS<T,std::enable_if_t<std::is_fundamental<T>::value>> {
        std::numeric_limits<T>::epsilon()
    };

template <typename T, typename Enabler = void> constexpr T huge;
template <typename T> constexpr T huge<T,std::enable_if_t<std::is_fundamental<T>::value>> {
        std::numeric_limits<T>::has_infinity
        ? std::numeric_limits<T>::infinity()
        : std::numeric_limits<T>::max()
    };

template <typename T, typename Enabler = void> constexpr T MAX;
template <typename T> constexpr T
    MAX<T,std::enable_if_t<std::is_fundamental<T>::value>> =
        std::numeric_limits<T>::max();
template <typename T> constexpr T
    MAX<T,std::enable_if_t<std::is_same<T,std::string>::value>> {""};

template <typename T, typename Enabler = void> constexpr T MIN;
template <typename T> constexpr T
    MIN<T,std::enable_if_t<std::is_fundamental<T>::value>> {
        std::numeric_limits<T>::min()
    };
template <typename T> constexpr T
    MIN<T,std::enable_if_t<std::is_same<T,std::string>::value>> {""};

template <typename T> constexpr T INF {std::numeric_limits<T>::infinity()};

// std arrays filled with common constants. ====================================

/// \brief produce std::array initialized to \p val.
/// \tparam T std::array template parameter.
/// \tparam N std::array template parameter.
/// \param val Value to which the array elements are set.
/// \return std::array initialized to \p val.
template <typename T, auto N> constexpr
auto filled_array( const T val )
{
    std::array<T,N> a {};
//    a.fill(val);        // needs c++20 to be constexpr
    for (auto& o : a)
        o = val;
    return a;
}

template <auto N> constexpr std::array<bool,N> falses {
    filled_array<bool,N>(false)
};
template <auto N> constexpr std::array<bool,N> trues {
    filled_array<bool,N>(true)
};
template <typename T, auto N> constexpr std::array<T,N> zeros {
    filled_array<T,N>(zero<T>)
};
template <typename T, auto N> constexpr std::array<T,N> ones {
    filled_array<T,N>(one<T>)
};
template <typename T, auto N> constexpr std::array<T,N> hundreds {
    filled_array<T,N>(static_cast<T>(100.L))
};
template <typename T, auto N> constexpr std::array<T,N> huges {
    filled_array<T,N>(huge<T>)
};
template <typename T, auto N> constexpr std::array<T,N> mhuges {
    filled_array<T,N>(-huge<T>)
};

// Range borders as two-element std arrays. ====================================

constexpr std::array<bool,2> bools {{false, true}};
template <typename T> constexpr std::array<T,2> zeroone {{zero<T>, one<T>}};
template <typename T> constexpr std::array<T,2> zerohuge {{zero<T>, huge<T>}};
template <typename T> constexpr std::array<T,2> onehuge {{one<T>, huge<T>}};
template <typename T> constexpr std::array<T,2> moneone {{-one<T>, one<T>}};
template <typename T> constexpr std::array<T,2> mhugehuge {{-huge<T>, huge<T>}};

// Range borders as two-element std vectors of arrays. =========================

template <auto N> const std::vector rangeBools {falses<N>, trues<N>};
template <auto N> const vecarr<bool,N> vecarrFT {falses<N>, trues<N>};
template <typename T, auto N> const vecarr<T,N> vecarr0H {zeros<T,N>, huges<T,N>};
template <typename T, auto N> const vecarr<T,N> vecarr01 {zeros<T,N>, ones<T,N>};
//template <typename T, auto N> constexpr std::vector<T,N> zeroshuges {zeros<real,N>, huges<real,N>};

// ANSI colors =================================================================

#define ANSI_RESET       "\x1B[0m"
#define ANSI_FG_BLACK    "\x1b[30m"
#define ANSI_FG_RED      "\x1B[31m"
#define ANSI_FG_GREEN    "\x1B[32m"
#define ANSI_FG_YELLOW   "\x1B[33m"
#define ANSI_FG_BLUE     "\x1B[34m"
#define ANSI_FG_MAGENTA  "\x1B[35m"
#define ANSI_FG_CYAN     "\x1B[36m"
#define ANSI_FG_WHITE    "\x1B[37m"
#define ANSI_BG_RED      "\x1b[41m"
#define ANSI_BG_GREEN    "\x1b[42m"
#define ANSI_BG_YELLOW   "\x1b[43m"
#define ANSI_BG_BLUE     "\x1b[44m"
#define ANSI_BG_MAGENTA  "\x1b[45m"
#define ANSI_BG_CYAN     "\x1b[46m"
#define ANSI_BG_WHITE    "\x1b[47m"
#define ANSI_BOLD_ON     "\x1b[1m"            
#define ANSI_BOLD_OFF    "\x1b[22m"        
#define ANSI_INVERSE_ON  "\x1b[7m"

}     // namespace Common
}    // namespace Utils


#endif // UTILS_COMMON_CONSTANTS_H
