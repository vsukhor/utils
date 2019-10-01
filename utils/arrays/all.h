/* ==============================================================================

 Copyright (C) 2009-2019, Valerii Sukhorukov, <vsukhorukov@yahoo.com>

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
* \file all.h
* \brief Generic array templates.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_ARRAYS_ARRAY_H
#define UTILS_ARRAYS_ARRAY_H

/// Library-wide.
namespace Utils {
/// Custom arrays.
namespace Arrays {

//===================================================================================
/// \brief A compile-time generator of std:iota-like array of consecutive integers [1, ..., N]
/// \details Adopted from https://stackoverflow.com/users/636019/ildjarn
/// at https://stackoverflow.com/questions/41660062/how-to-construct-an-stdarray-with-index-sequence
namespace iota_array {
  template<typename T, T... Ns>
  constexpr std::array<T, sizeof...(Ns)> make_iota_array(T const offset,
  														 std::integer_sequence<T, Ns...>) noexcept
  {
    return {{(Ns + offset)...}};
  }
}

template<typename T, T N>
constexpr auto make_iota_array(T const offset = {}) noexcept
{
  static_assert(N >= T{}, "no negative sizes");
  return iota_array::make_iota_array<T>(offset, std::make_integer_sequence<T, N>{});
}

//===================================================================================

/// \brief Generic array template.
/// \tparam N Array length.
/// \tparam T Type of the elements.
/// \tparam Enabler SFINAE Enabler
template <unsigned N, typename T, typename Enabler=void>
class array {};

}	// namespace Arrays
}	// namespace Utils

#include "array2.h"
#include "array3.h"
#include "array4.h"
#include "arrayN.h"

namespace Utils {
namespace Arrays {

// Type abbreviations.
template<typename T> using A2 = Arrays::array<2,T>;
template<typename T> using A3 = Arrays::array<3,T>;
template<typename T> using A4 = Arrays::array<4,T>;

}	// namespace Arrays
}	// namespace Utils

#endif // UTILS_ARRAYS_ARRAY_H
