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

============================================================================== */

/**
 * \file exceptions.h
 * \brief std::exception-derived classes for checking confuguration parameters.
 * \details Implementation of std::exception-derived classes for checking
 * the parameters read from confuguration files.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_EXCEPTIONS_H
#define UTILS_CONFIG_EXCEPTIONS_H

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "msgr.h"

/// Exeption and exit handlers.
namespace utils::common::exceptions {

/**
 * \brief Simple process termination as a global function.
 * \details Outputs message \p msg to \a Msgr.
 * \param msg Message to output.
 * \see Msgr
*/
int simple(const std::string& msg,
           Msgr* msgr);


// Simple exceptions ===========================================================

/**
* \brief Uncomplicated process termination with exception.
*/
class Simple
    : public std::exception {

public:

    /// \brief Default constructor.
    Simple() = default;

    /**
    * \brief Constructor for printing to a log record.
    * \details Outputs message \p msg to \a Msgr for standard
    * output (if \p msgr is nullptr).
    * \param msg Message to output.
    * \param msgr Output message processor.
    */
    explicit Simple(const std::string& msg,
                    Msgr* msgr);
};

}    // namespace utils::common::exceptions

#endif // UTILS_CONFIG_EXCEPTIONS_H
