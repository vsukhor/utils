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
* \file msgr.cpp
* Implements class Msgr.
* \author Valerii Sukhorukov
*/

#include "msgr.h"

/// Library-wide.
namespace Utils {
/// General stuff.
namespace Common {


Msgr::
Msgr( outstream* so,
      logstream* sl,
      const int precision
    )
    : so {so}
    , sl {sl}
{
    set_formats(precision);
}


void Msgr::
set_formats( const int precision ) noexcept
{
    if (so) {
        so->precision(precision);
        so->setf(std::ios::scientific);
    }
    if (sl) {
        sl->precision(precision);
        sl->setf(std::ios::scientific);
    }
}


void Msgr::
exit( const std::string& s ) const noexcept
{    
    print<true>(s);
    ::exit(EXIT_FAILURE);
}


void Msgr::
exit( const char *fmt, ... ) noexcept
{    
    va_list va;
    va_start(va, fmt);
    const auto n = vsprintf(buf, fmt, va);
    va_end(va);
    const auto s = std::string(buf).substr(0, std::size_t(n));
    if (sl) prn(sl, s, 1);
    if (so) prn(so, s, 1);

    ::exit(EXIT_FAILURE);
}

}    // namespace Common
}    // namespace Utils
