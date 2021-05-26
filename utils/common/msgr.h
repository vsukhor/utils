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
* \file msgr.h
* \brief Contains class Msgr.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_COMMON_MSGR_H
#define UTILS_COMMON_MSGR_H

#include <array>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <vector>

/// General stuff.
namespace Utils::Common {

/**
* \class Msgr msgr.h
* \brief Convenient formatted text output to the screen and a logfile.
* \details Implements convenience class for formatted text output
* to std::cout and to a logfile.
*/
class Msgr {

    using outstream = std::ostream;
    using logstream = std::ofstream;


public:

    outstream* so {};  ///< Screen out stream.
    logstream* sl {};  ///< Logfile stream.

    /// Default constructor.
    Msgr() = default;

    /**
    * \brief Constructor.
    * \param so Screen out stream.
    * \param sl File out stream.
    * \param precision Precision of real numbers.
    */
    explicit Msgr(
        outstream* so,
        logstream* sl,
        int precision
    );
    
    /**
    * \brief Set formatting parameters.
    * \param precision Precision of real numbers
    */
    void set_formats(int precision) noexcept;


    /**
    * \brief Print std::array.
    * \tparam V Data type of array elements.
    * \tparam N Number of array elements.
    * \param name Name/title.
    * \param v Array data.
    */
    template <typename V, auto N>
    void print_array( const std::string& name,
                      const std::array<V,N>& v
                    ) const noexcept;

    /**
    * \brief Print named std::vector .
    * \tparam V Data type of vector elements.
    * \param name Name/title.
    * \param v Vector data.
    */
    template <typename V>
    void print_vector( const std::string& name,
                       const std::vector<V>& v
                     ) noexcept;

    /**
    * \brief Print std::string.
    * \tparam endline Finish with line end.
    * \param s String to print.
    */
    template <bool endline=true>
    void print(const std::string& s) const noexcept;

    /**
    * \brief Print std::string out and exit.
    * \param s String to print.
    */
    void exit(const std::string& s) const noexcept;

private:

    /**
    * \brief Check that the stream used is valid.
    * \tparam S Stream type.
    */
    template <typename S>
    static constexpr auto is_valid_stream() noexcept;

    /**
    * \brief Print to a stream \p io.
    * \tparam IO Stream type.
    * \param v String to print.
    * \param endline Specifies if the line end should be added.
    */
    template <typename IO>
    void prn(IO* io,
             const std::string& v,
             bool endline
             ) const noexcept;
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S>
constexpr auto Msgr::
is_valid_stream() noexcept
{
    return std::is_same_v<S, outstream> ||
           std::is_same_v<S, logstream>;
}


template <typename IO> inline
void Msgr::
prn(
    IO* io,
    const std::string& v,
    const bool endline
) const noexcept
{ 
    static_assert(is_valid_stream<IO>(),
                  "Stream type used in Msgr is not valid");

    *io << v << " ";
    if (endline) *io << std::endl;
}


template <bool endline> inline
void Msgr::
print( const std::string& s ) const noexcept
{    
    if (sl) prn(sl, s, endline);
    if (so) prn(so, s, endline);
}


template <typename V, auto N>
void Msgr::
print_array(
    const std::string& name,
    const std::array<V,N>& v
) const noexcept
{
    print<false>(name+"[]:  ");
    for (const auto o : v)
        print<false>(std::to_string(o));
    print<true>("");
}


template <typename V>
void Msgr::
print_vector(
    const std::string& name,
    const std::vector<V>& v
) noexcept
{
    print<false>(name+"[]:  ");
    for (const auto o : v)
        print<false>(std::to_string(o));
}

}  // namespace Utils::Common

#endif // UTILS_COMMON_MSGR_H
