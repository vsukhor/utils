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
 * \file exceptions.h
 * \brief std::exception-derived classes for checking the parameters read from confuguration files.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_EXCEPTIONS_H
#define UTILS_CONFIG_EXCEPTIONS_H

#include <iostream>
#include <string>
#include <algorithm>
#include <exception>
#include <vector>

#include "msgr.h"

namespace Utils {
namespace Common {
namespace Exceptions {

/**
* Simple process termination as a global function.
* Outputs message \p msg to \a Msgr.
* \param msg Message to output.
* \return EXIT_FAILURE
* \see Msgr
*/
int simple(const std::string& msg,
		   Msgr* msgr=nullptr);

// Simple Exceptions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
* Simple process termination with exception.
*/
class Simple
	: public std::exception {

public:

	/**
	* Default constructor.
	*/
	Simple() = default;

	/**
	* Constructor for printing to a log record.
	* Outputs message \p msg to \a Msgr sor tandard output (if \p msgr is nullptr).
	* \param msg Message to output.
	* \see Msgr
	*/
	explicit Simple(const std::string& msg,
					Msgr* msgr=nullptr);

	/**
	* Constructor template for message with parameter \u.
	* Outputs message \p msg to \a Msgr.
	* \param msg Message to output.
	* \param u Message parameter.
	* \tparam T Type of the essage parameter.
	* \see Msgr
	*/
	template<typename T>
	explicit Simple(const std::string& msg,
					const T& u,
					Msgr* msgr=nullptr);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename T>
Simple::
Simple(
		const std::string& msg,
		const T& u,
		Msgr* msgr )
	: std::exception {}
{
	if (msgr)
		msgr->print(msg, u);
	else
		std::cout << msg << " " << u << std::endl;
}


}	// namespace Exceptions
}	// namespace Common
}	// namespace Utils

#endif // UTILS_CONFIG_EXCEPTIONS_H
