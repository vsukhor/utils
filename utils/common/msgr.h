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
* \file msgr.h
* Contains class Msgr.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_COMMON_OEL_H
#define UTILS_COMMON_OEL_H

#include <type_traits>
#include <array>
#include <fstream>
#include <iostream>
#include <stdarg.h>

namespace Utils {
namespace Common {

/**
* \class Gillespie gillespie.h
* \brief Convenient formatted text output to the screen and a logfile.
* Implements convenience class for formatted text output to std::cout and to a logfile.
*/
class Msgr {

	using outstream = std::ostream;
	using logstream = std::ofstream;


public:

	outstream* so {};		///< Screen out stream.
	logstream* sl {};		///< Logfile stream.

	/**
	* Default constructor.
	*/
	Msgr() = default;

	/**
	* Constructor.
	* \param so Screen out stream.
	* \param sl File out stream.
	* \param precision Precision of real numbers.
	*/
	explicit Msgr(outstream* so,
				  logstream* sl,
		  		  const int precision=6 );
	
	/**
	* Set formatting parameters.
	* \param precision Precision of real numbers
	*/
	void set_formats(const int precision) noexcept;


	/**
	* Print std::array out.
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
	* Print std::string out.
	* \tparam endline Finish with line end.
	* \param s String to print.
	*/
	template <bool endline=true>
	void print(const std::string& s) const noexcept;

	/**
	* Print to formatted out (variadic).
	* \tparam endline Finish with line end.
	* \param fmt Formatting.
	*/
	template <bool endline=true>
	void print(const char *fmt, ...) noexcept;
		
	/**
	* Print std::string out and exit.
	* \param s String to print.
	*/
	void exit(const std::string& s) const noexcept;

	/**
	* Print to formatted out (variadic) and exit.
	* \param fmt Formatting.
	*/
	void exit(const char *fmt, ...) noexcept;
	
private:

	char buf [4096];	///< Buffer.

	/**
	* Check that the stream used is valid.
	* \tparam S Stream type.
	*/
	template <typename S>
	static constexpr auto is_valid_stream() noexcept;

	/**
	* Print to a stream \p io.
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

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename S>
constexpr auto Msgr::
is_valid_stream() noexcept {
	return std::is_same_v<S, outstream> ||
		   std::is_same_v<S, logstream>;
}

template <typename IO> inline
void Msgr::
prn( IO* io,
	const std::string& v,
	bool endline ) const noexcept
{ 
	static_assert(is_valid_stream<IO>(), "Stream type used in Msgr is not valid");

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

template <bool endline>
void Msgr::
print( const char *fmt, ... ) noexcept
{
	va_list va;
	va_start(va, fmt);
	const auto n = vsprintf(buf, fmt, va);
	va_end(va);
	const auto s = std::string(buf).substr(0, static_cast<unsigned long>(n));
	if (sl) prn(sl, s, endline);
	if (so) prn(so, s, endline);
}

template <typename V, auto N>
void Msgr::
print_array( const std::string& name,
		  	 const std::array<V,N>& v ) const noexcept
{
	print<false>(name+"[]:  ");
	for (const auto o : v)
		print<false>(std::to_string(o));
	print("");
}

}	// namespace Common
}	// namespace Utils

#endif // UTILS_COMMON_OEL_H
