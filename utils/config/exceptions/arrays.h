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
 * \file arrays.h
 * \brief Exception handler for parameters of std::array type.
 * \details A class template specialization for checking parameters
 * of type std::array read from confuguration files.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_EXCEPTIONS_ARRAYS_H
#define UTILS_CONFIG_EXCEPTIONS_ARRAYS_H

#include <algorithm>
#include <array>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

#include "../../msgr.h"
//#include "../parameter/base.h"
#include "generic.h"

namespace utils::config::exceptions {

/**
 * \brief 'Parameter out of range' exception for continuous arithmetic arrays.
 * \details Partial template specialization for for 'Parameter out of range'
 * exception involving continuous arithmetic std arrays.
 * \tparam T Parameter type: must be std::is_arithmetic.
 * \tparam W Length of the std::array.
 */
template<typename T,
         std::size_t W>
class ParOutOfRange<std::array<T, W>, false,
                    std::enable_if_t<std::is_arithmetic_v<T>>>
    : public std::exception {

    using Q = std::array<T, W>;

public:

    /**
     * \brief Constructor.
     * \details Upon throwing, prints message generated by \a generate_message
     * to \p msgr.
     * \param name Name of the parameter that triggered the exception.
     * \param p Parameter value (out of range).
     * \param r Acceptable range of parameter values.
     * \param msgr Pointer to \a Msgr used for the output.
     * \see Msgr
     */
    explicit ParOutOfRange(
        const std::string& name,
        const Q& p,
        const std::vector<Q>& r,
        Msgr* msgr=nullptr
    );
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename T,
         std::size_t W>
ParOutOfRange<std::array<T, W>, false,
              std::enable_if_t<std::is_arithmetic_v<T>>>::
ParOutOfRange(
    const std::string& name,
    const Q& p,
    const std::vector<Q>& r,
    Msgr* msgr
)
    : std::exception {}
{
    ASSERT((r.size() == 2 && std::all_of(r.begin(), r.end(),
            [&p](const Q& a){return a.size()==p.size();})),
           "Incorrect r size in ParOutOfRangeException");

    auto print = [](const Q& a) {
        std::string w{"{ "};
        for (const auto& o : a)
            w += "[" + std::to_string(o[0]) + " " + std::to_string(o[1]) + "] ";
        return w+"}";
    };

    const std::string s =
        "Error in conf specification for parameter '" +
        name + "' = " + print(p) + " :" + "\n\tthe value " +
        " is outside the acceptable range " + print(r);

    if (msgr != nullptr)
        msgr->print<true>(s);
    else
        std::cerr << s << std::endl;
}

}  // namespace utils::config::exceptions

#endif // UTILS_CONFIG_EXCEPTIONS_ARRAYS_H
