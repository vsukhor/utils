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
 * \file scalars.h
 * \brief Exception handler for parameters of scalar type.
 * \details A class template specialization for checking parameters of scalar
 * arithmetic type read from confuguration files.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_EXCEPTIONS_SCALARS_H
#define UTILS_CONFIG_EXCEPTIONS_SCALARS_H

#include <algorithm>
#include <exception>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "../../common/msgr.h"
#include "../parameter/base.h"

/// Exceptions namespace.
namespace utils::config::exceptions {

/**
* \brief 'Parameter out of range' exception involving discrete scalars.
* \details Partial template specialization for for 'Parameter out of range'
* exception involving discrete scalars.
* \tparam Q Parameter type: must be std::is_arithmetic.
*/
template <typename Q>
class ParOutOfRange<Q, true,        // discrete case
                    std::enable_if_t<std::is_arithmetic_v<Q>>>
    : public std::exception {

public:

    /**
    * Constructor.
    * \details Upon throwing, prints message generated by \a generate_message
    * to \p msgr or standard output (if \p msgr is \a nullptr).
    * \param name Name of the parameter that triggered the exception.
    * \param p Parameter value (out of range).
    * \param r Acceptable range of parameter values.
    * \see Msgr
    */
    explicit ParOutOfRange(
            const std::string& name,
            const Q& p,
            const std::vector<Q>& r,
            Msgr* msgr=nullptr);

private:

    const std::string message;  ///< message to print

    /**
    * Message-generating function.
    * \details Outputs message \p msg to \a Msgr.
    * \param name Name of the parameter that triggered the exception.
    * \param p Parameter value (out of range).
    * \param r Acceptable range of parameter values.
    */
    std::string generate_message(
            const std::string& name,
            const Q& p,
            const std::vector<Q>& r);
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Q>
ParOutOfRange<Q, true,        // discrete case
              std::enable_if_t<std::is_arithmetic_v<Q>>>::
ParOutOfRange(
        const std::string& name,
        const Q& p,
        const std::vector<Q>& r,
        Msgr* msgr)
    : std::exception {}
    , message {generate_message(name, p, r)}
{
    if (msgr != nullptr)
        msgr->print(message);
    else
        std::cout << message;
        
    throw this;
}

template <typename Q>
std::string ParOutOfRange<Q, true,        // discrete case
            std::enable_if_t<std::is_arithmetic_v<Q>>>::
generate_message(const std::string& name,
                 const Q& p,
                 const std::vector<Q>& r)
{
    XASSERT(r.size() == 2, "Incorrect r size in ParOutOfRangeException");

    auto print = [](const std::vector<Q>& a) {
        std::string w {"{ "};
        for(const auto o : a)
            w += std::to_string(o)+" ";
        return w+"}";
    };

    return "Error in conf specification for parameter '"+name+
           "' = "+std::to_string(p)+" :"+"\n\tthe value provided "+
           " is outside thw acceptable range "+
           print(r);
}


////////////////////////////////////////////////////////////////////////////////
/**
* \brief 'Parameter out of range' exception.
* \details Partial template specialization for for 'Parameter out of range'
* exception involving continuous arithmetic scalars.
* \tparam Q Parameter type: must be std::is_arithmetic.
*/
template <typename Q>
class ParOutOfRange<Q, false,        // continuous case
                    std::enable_if_t<std::is_arithmetic_v<Q>>>
    : public std::exception {

public:

    /**
    * \brief Constructor.
    * \details Upon throwing, prints message generated by \a generate_message
    * to \p msgr.
    * \param name Name of the parameter that triggered the exception.
    * \param p Parameter value (out of range).
    * \param r Acceptable range of parameter values.
    * \param msgr Pointer to \a Msgr used for the output.
    */
    explicit ParOutOfRange(
            const std::string& name,
            const Q& p,
            const std::vector<Q>& r,
            Msgr* msgr=nullptr);

    /**
    * \brief Constructor.
    * \details Upon throwing, prints message generated by \a generate_message
    * to \p msgr.
    * \param name Name of the parameter that triggered the exception.
    * \param p Parameter value (out of range).
    * \param r Acceptable range of parameter values.
    * \param msgr Pointer to \a Msgr used for the output.
    */
    explicit ParOutOfRange(
            const std::string& name,
            const Q& p,
            const std::array<Q,2>& r,
            Msgr* msgr=nullptr);

private:

    const std::string message;  ///< Message to print.

    /**
    * \brief Message-generating function.
    * \details Outputs message \p msg to \a Msgr.
    * \param name Name of the parameter that triggered the exception.
    * \param p Parameter value (out of range).
    * \param r Acceptable range of parameter values.
    * \tparam R Type of the range specificator.
    */
    template <typename R>
    std::string generate_message(const std::string& name,
                                 const Q& p,
                                 const R& r);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Q>
ParOutOfRange<Q, false,        // continuous case
              std::enable_if_t<std::is_arithmetic_v<Q>>>::
ParOutOfRange(
        const std::string& name,
        const Q& p,
        const std::vector<Q>& r,
        Msgr* msgr)
    : std::exception {}
    , message {generate_message(name, p, r)}
{
    if (msgr != nullptr)
        msgr->print(message);
    else
        std::cerr << message;
        
    throw *this;
}

template <typename Q>
ParOutOfRange<Q, false,        // continuous case
                std::enable_if_t<std::is_arithmetic_v<Q>>>::
ParOutOfRange(
        const std::string& name,
        const Q& p,
        const std::array<Q,2>& r,
        Msgr* msgr)
    : std::exception {}
    , message {generate_message(name, p, r)}
{
    if (msgr != nullptr)
        msgr->print(message);
    else
        std::cerr << message;

    throw *this;
}

template <typename Q>
template <typename R>
std::string ParOutOfRange<Q, false,        // continuous case
                          std::enable_if_t<std::is_arithmetic_v<Q>>>::
generate_message(const std::string& name,
                 const Q& p,
                 const R& r)
{
    XASSERT(r.size() == 2, "Incorrect r size in ParOutOfRangeException");
    return "Error in conf specification for parameter '"+name+
           "' = "+std::to_string(p)+" :"+"\n\tthe value "+
           " is outside the acceptable range "+
           "[ "+std::to_string(r[0])+", "+std::to_string(r[1])+" ].\n";
}

}   // namespace utils::config::exceptions

#endif // UTILS_CONFIG_EXCEPTIONS_SCALARS_H
