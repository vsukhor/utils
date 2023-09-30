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

/**
 * \file base.h
 * \brief Contains the base class for configuration parameters.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_EXCEPTIONS_GENERIC_H
#define UTILS_CONFIG_EXCEPTIONS_GENERIC_H

#include <exception>

/// Processing of configuration files.
namespace utils::config::exceptions {

/**
 * \brief Generic template for 'Parameter out of range' exception.
 * \tparam Q Parameter type.
 * \tparam isDiscrete Specifies if the parameter takes discrete values only.
 * \tparam Enabler SFINAE enabler.
 */
template<typename Q, 
         bool isDiscrete, 
         typename Enabler = void>
class ParOutOfRange
    : public std::exception
{};

}  // namespace exceptions

#endif  // UTILS_CONFIG_EXCEPTIONS_GENERIC_H