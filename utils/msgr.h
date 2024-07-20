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
 
================================================================================
*/

/**
 * \file msgr.h
 * \brief Contains class Msgr.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_MSGR_H
#define UTILS_MSGR_H

#include <array>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdarg.h>
#include <string>
#include <type_traits>
#include <utility>   // forward
#include <vector>

/// Library outer namespace.
namespace utils {


struct Colorcodes {

static constexpr const char* BOLD        {"\033[1m"};
static constexpr const char* RESET       {"\033[0m"};
static constexpr const char* BLACK       {"\033[30m"};          // Black
static constexpr const char* RED         {"\033[31m"};          // Red 
static constexpr const char* GREEN       {"\033[32m"};          // Green 
static constexpr const char* YELLOW      {"\033[33m"};          // Yellow 
static constexpr const char* BLUE        {"\033[34m"};          // Blue 
static constexpr const char* MAGENTA     {"\033[35m"};          // Magenta 
static constexpr const char* CYAN        {"\033[36m"};          // Cyan 
static constexpr const char* WHITE       {"\033[37m"};          // White 
static constexpr const char* BOLDBLACK   {"\033[1m\033[30m"};   // Bold Black 
static constexpr const char* BOLDRED     {"\033[1m\033[31m"};   // Bold Red 
static constexpr const char* BOLDGREEN   {"\033[1m\033[32m"};   // Bold Green 
static constexpr const char* BOLDYELLOW  {"\033[1m\033[33m"};   // Bold Yellow 
static constexpr const char* BOLDBLUE    {"\033[1m\033[34m"};   // Bold Blue 
static constexpr const char* BOLDMAGENTA {"\033[1m\033[35m"};   // Bold Magenta 
static constexpr const char* BOLDCYAN    {"\033[1m\033[36m"};   // Bold Cyan 
static constexpr const char* BOLDWHITE   {"\033[1m\033[37m"};   // Bold White 

const std::string sBOLD        {BOLD};
const std::string sRESET       {RESET};
const std::string sBLACK       {BLACK};          
const std::string sRED         {RED};         
const std::string sGREEN       {GREEN};          
const std::string sYELLOW      {YELLOW};         
const std::string sBLUE        {BLUE};          
const std::string sMAGENTA     {MAGENTA};         
const std::string sCYAN        {CYAN};          
const std::string sWHITE       {WHITE};          
const std::string sBOLDBLACK   {BOLDBLACK};   
const std::string sBOLDRED     {BOLDRED};   
const std::string sBOLDGREEN   {BOLDGREEN};   
const std::string sBOLDYELLOW  {BOLDYELLOW};   
const std::string sBOLDBLUE    {BOLDBLUE};   
const std::string sBOLDMAGENTA {BOLDMAGENTA};   
const std::string sBOLDCYAN    {BOLDCYAN};   
const std::string sBOLDWHITE   {BOLDWHITE};  


static constexpr std::array all {
    BOLD, RESET, 
    BLACK,     RED,     GREEN,     YELLOW,     BLUE,     MAGENTA,     CYAN,     WHITE,
    BOLDBLACK, BOLDRED, BOLDGREEN, BOLDYELLOW, BOLDBLUE, BOLDMAGENTA, BOLDCYAN, BOLDWHITE
};

const std::array<std::string, all.size()> sall {
    sBOLD, sRESET, 
    sBLACK,     sRED,     sGREEN,     sYELLOW,     sBLUE,     sMAGENTA,     sCYAN,     sWHITE,
    sBOLDBLACK, sBOLDRED, sBOLDGREEN, sBOLDYELLOW, sBOLDBLUE, sBOLDMAGENTA, sBOLDCYAN, sBOLDWHITE
};

};  // struct Colorcodes

constexpr Colorcodes colorcodes;


// /////////////////////////////////////////////////////////////////////////////

/**
 * \class Msgr msgr.h
 * \brief Convenient formatted text output to the screen and a logfile.
 * \details Implements convenience class for formatted text output
 * to std::cout and to a logfile.
 */
class Msgr {

public:

    using outstream = std::ostream;
    using logstream = std::ofstream;


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
    template<typename V, auto N>
    void print_array(
        const std::string& name,
        const std::array<V,N>& v
    ) const noexcept;

    /**
     * \brief Print named std::vector .
     * \tparam V Data type of vector elements.
     * \param name Name/title.
     * \param v Vector data.
     */
    template<typename V>
    void print_vector(
        const std::string& name,
        const std::vector<V>& v
    ) const noexcept;

    /**
     * \brief Print an data series of of various types.
     * \details Based on proposition stackoverflow.com proposition by geza:
     * https://stackoverflow.com/questions/51647834/printf-like-utility-in-c-without-format-specifier/51648068#51648068
     * \tparam endline Finish with line end.
     * \param values Values to print.
     */
    template<bool endline=true, typename... T>
    void print(T... values) const;

    /**
     * \brief Print an data series of of various types and exit the process.
     * \details Based on proposition stackoverflow.com proposition by geza:
     * https://stackoverflow.com/questions/51647834/printf-like-utility-in-c-without-format-specifier/51648068#51648068
     * \param values Values to print.
     */
    template<typename... T>
    void exit(T... values);

    template<typename... T>
    Msgr& operator<<(T... values);
    
private:

    /**
     * \brief Check that the stream used is valid.
     * \tparam S Stream type.
     */
    template<typename S>
    static constexpr auto is_valid_stream() noexcept;

    /**
     * \brief Removes color codes from the stream directed to a file.
     */
    auto rm_colorcode(std::string&& v) const noexcept -> std::string;

    /**
     * \brief Print to a stream \p io.
     * \tparam IO Stream type.
     * \param v String to print.
     * \param endline Specifies if the line end should be added.
     */
    template<typename IO>
    void prn(
        IO* io,
        std::string&& v,
        bool endline
    ) const noexcept;
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

inline
Msgr::
Msgr(
    outstream* so,
    logstream* sl,
    int precision
)
    : so {so}
    , sl {sl}
{
    set_formats(precision);
}


inline
void Msgr::
set_formats(const int precision) noexcept
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


template<typename S>
constexpr auto Msgr::
is_valid_stream() noexcept
{
    return std::is_same_v<S, outstream> ||
           std::is_same_v<S, logstream>;
}


inline
auto Msgr::
rm_colorcode(std::string&& s) const noexcept -> std::string
{
    for (const auto& p: colorcodes.sall) {
        auto n = p.length(); 
        for (auto i = s.find(p); i != std::string::npos; i = s.find(p)) 
            s.erase(i, n);
    }

    return s;
}


template<typename IO>
void Msgr::
prn(
    IO* io,
    std::string&& v,
    const bool endline
) const noexcept
{ 
    static_assert(is_valid_stream<IO>(),
                  "Stream type used in Msgr is not valid");

    if constexpr (std::is_same_v<IO, logstream>) 
        v = rm_colorcode(std::forward<std::string>(v));
    
    *io << v << " ";
    if (endline)
        *io << std::endl;
    io->flush();
}


template<typename V, auto N>
void Msgr::
print_array(
    const std::string& name,
    const std::array<V,N>& v
) const noexcept
{
    print<false>(name+"[]:  ");
    for (const auto o : v)
        print<false>(o);
    print("");
}


template<typename V>
void Msgr::
print_vector(
    const std::string& name,
    const std::vector<V>& v
) const noexcept
{
    print<false>(name+"[]:  ");
    for (const auto o : v)
        print<false,V>(o);
    print<true>("");
}


template<bool endline,
         typename... T>
void Msgr::
print(T... values) const
{
    std::ostringstream s;

    (s << ... << values);
    s.flush();

	if (sl) prn(sl, s.str(), endline);
	if (so) prn(so, s.str(), endline);
}


template<typename... T>
void Msgr::
exit(T... values)
{
    std::ostringstream s;

    (s << ... << values);
    s.flush();

	if (sl) prn(sl, s.str(), true);
	if (so) prn(so, s.str(), true);

    ::exit(EXIT_FAILURE);
}


template<typename... T>
Msgr& Msgr::
operator<<(T... values)
{
    print<false>(values...);
    
    return *this;
}


}  // namespace utils

#endif // UTILS_MSGR_H
