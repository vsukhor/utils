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
 * \file exceptions.cpp
 * \brief Implementation of std::exception-derived classes for checking the parameters read from confuguration files.
 + \author Valerii Sukhorukov
 */

#include "exceptions.h"

/// Library-wide.
namespace Utils {
/// General stuff.
namespace Common {
/// \brief Exeption and exit handlers.
namespace Exceptions {

/// \brief Print a message and exit.
/// \param msg Message text.
/// \param msg The reporter.
int simple( const std::string& msg,
			Msgr* msgr )
{
	if (msgr)
		msgr->print<true>(msg);
	else
		std::cout << msg << std::endl;

	return EXIT_FAILURE;
}

// Simple xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
Simple::
Simple(
		const std::string& msg,
		Msgr* msgr )
	: std::exception {}
{
	if (msgr)
		msgr->print<true>(msg);
	else
		std::cout << msg << std::endl;
}

}	// namespace Exceptions
}	// namespace Common
}	// namespace Utils
