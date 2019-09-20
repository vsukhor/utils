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

#ifndef UTILS_COMMON_ROOT_FINDER_H
#define UTILS_COMMON_ROOT_FINDER_H

#include "msgr.h"
#include "misc.h"

namespace Utils {
namespace Common {

template <int I, typename T>
class RootFinder {

public:

	/**
	* \brief Constructor.
	* \param msgr Output printing.
	*/
	explicit RootFinder(Msgr& msgr);

	 T rtflsp_vTvTvTT( T (*func)(const std::array<T,I>&, const std::array<T,I>&, const std::array<T,I>&, const T),
					  const T x1, const T x2, const T xacc,
					  const std::array<T,I>& arg1, const std::array<T,I>& arg2, const std::array<T,I>& arg3 ) const;

	T rtflsp_monotGr_TvTvTT( T (*func)( const std::array<T,I>&, const std::array<T,I>&, const std::array<T,I>&, const T ),
							const T xl, const T xh, const T xacc,
							const std::array<T,I>& arg1, const std::array<T,I>& arg2, const std::array<T,I>& arg3 ) const;
private:

	Msgr& 				  msgr;					///< Output handling facility.
	static constexpr uint MAXIT {100'000'000};	///< Max. number of iterations.
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

RootFinder( Msgr& msgr )
	: msgr {msgr}
{}

template <int I, typename T>
T RootFinder<I,T>::
rtflsp_vTvTvTT( T (*func)(const std::array<T,I>&, const std::array<T,I>&, const std::array<T,I>&, const T),
				  const T x1, const T x2, const T xacc,
				  const std::array<T,I>& arg1, const std::array<T,I>& arg2, const std::array<T,I>& arg3 ) const
{
	// Using the false position method, find the root of a function func known to lie between x1 and x2.
	// The root, returned as rtflsp, is refined until its accuracy is ±xacc.
	// void nrerror(char error_text[]);

	T fl = (*func)( arg1, arg2, arg3, x1 );
	T fh = (*func)( arg1, arg2, arg3, x2 );			// Be sure the interval brackets a root. 
	if (fl*fh > zero<T>)
		msgr.exit("Error: Root must be bracketed in rtflsp");

	T xl, xh;
	if (fl < zero<T>) {
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2; 
		xh = x1; 
		const auto swap {fl};
		fl = fh;
		fh = swap; 
	}
	auto dx {xh - xl};
	T del {}, rtf {};
	for (uint j=1; j<=MAXIT; j++) {				// False position loop. Increment with respect to latest value.
		rtf = xl + dx * fl / (fl - fh);
		const auto f = (*func)(arg1, arg2, arg3, rtf);
		if (f < zero<T>) {						// Replace appropriate limit.
			del = xl - rtf; 
			xl = rtf; 
			fl = f;
		} else { 
			del = xh - rtf;
			xh = rtf;
			fh = f; 
		}
		dx = xh - xl;
		if (std::abs(del) < xacc || f == zero<T>)
			return rtf;		// Convergence.
	}
	msgr.exit("Maximum number of iterations exceeded in rtflsp: del "+STR(del)+", f "+STR((*func)(arg1, arg2, arg3, rtf))); // Never get here.
	return huge<T>;
}

template <int I, typename T>
T RootFinder<I,T>::
rtflsp_monotGr_TvTvTT( T (*func)( const std::array<T,I>&, const std::array<T,I>&, const std::array<T,I>&, const T ),
						const T xl, const T xh, const T xacc,
						const std::array<T,I>& arg1, const std::array<T,I>& arg2, const std::array<T,I>& arg3 ) const
{
	// Using the false position method, find the root of a function func known to lie between x1 and x2.
	// The root, returned as rtflsp, is refined until its accuracy is ±xacc.
	// void nrerror(char error_text[]);

	auto fl = (*func)(arg1, arg2, arg3, xl);
	while (fl >= zero<T>) {
		xl /= three<T>;
		fl = (*func)(arg1, arg2, arg3, xl);
	}

	T fh = (*func)(arg1, arg2, arg3, xh);			// Be sure the interval brackets a root.
	while (fh <= zero<T>) {
		xh /= three<T>;
		fh = (*func)(arg1, arg2, arg3, xh);
	}

	auto dx {xh - xl};
	T del {}, rtf {};
	for (uint j=1; j<=MAXIT; j++) {				// False position loop. Increment with respect to latest value.
		rtf = xl + dx * fl / (fl - fh);
		const auto f = (*func)(arg1, arg2, arg3, rtf);
		if (f < zero<T>) {						// Replace appropriate limit.
			del = xl - rtf; 
			xl = rtf; 
			fl = f;
		} else { 
			del = xh - rtf;
			xh = rtf;
			fh = f; 
		}
		dx = xh - xl;
		if (std::abs(del) < xacc || f == zero<T>)
			return rtf;		// Convergence.
	}
	return msgr.exit("Maximum number of iterations exceeded in rtflsp: del "+STR(del)+", f "+STR((*func)(arg1, arg2, arg3, rtf))); // Never get here.
}

}	// namespace Common
}	// namespace Utils

#endif	// UTILS_COMMON_ROOT_FINDER_H















