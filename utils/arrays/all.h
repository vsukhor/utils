/* =============================================================================

 Copyright (C) 2009-2023 Valerii Sukhorukov. All Rights Reserved.

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

/***
 * \file all.h
 * \brief Generic array templates.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_ARRAYS_ARRAY_H
#define UTILS_ARRAYS_ARRAY_H

/// Custom arrays.
#include "../constants.h"
#include "array2.h"
#include "array3.h"
#include "array4.h"
#include "arrayN.h"

namespace utils::arrays {


// Type abbreviations.
template<arithmetic T> using A2 = arrays::array<2, T>;
template<arithmetic T> using A3 = arrays::array<3, T>;
template<arithmetic T> using A4 = arrays::array<4, T>;


}  // namespace utils::arrays

#endif // UTILS_ARRAYS_ARRAY_H
