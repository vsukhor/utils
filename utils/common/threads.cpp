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
* \file threads.h
* Implements class Threads.
* \author Valerii Sukhorukov
*/

#include "threads.h"

namespace Utils {
namespace Common {

Threads::
Threads(const szt offset,
		const szt size,
		const ulong omittedBoundaries,
		const Weights wht,
		const ulong nThreads )
	: num {nThreads}
{
	szt threadsSupported {std::thread::hardware_concurrency()};
	if (nThreads > threadsSupported) {
		std::cout << "Warning: thread number set exceeds the CPU concurrency: " << nThreads << " " << threadsSupported << std::endl; 
		exit(0); 
	}
	if (nThreads < 1) {
		std::cout << "Error in Threads: nThreads provided is not supprted " << nThreads << std::endl;
		exit(0); 
	}
	szt w {size - 2*omittedBoundaries};
	const szt rest {w % num};
	w -= rest;
	
	chunkSize.resize(num);
	if (     wht == Weights::Equal)		   set_chunks_equal(w, rest);
	else if (wht == Weights::CircleCenter) set_chunks_circular(w, rest);	// arc sector area: A = r*r*phi/2
	else if (wht == Weights::TriangleDecr) set_chunks_triangleDecr(size - 2*omittedBoundaries);	
	else { 
		std::cout << "Error in Threads: Weight type is not defined" << std::endl; 
		exit(0); 
	}
	i1.resize(num);
	i2.resize(num);
	i1[0] = offset + omittedBoundaries;
	i2[0] = i1[0] + chunkSize[0];
	for (szt ith=1; ith<num; ith++) {
		i1[ith] = i1[ith-1] + chunkSize[ith-1];
		i2[ith] = i2[ith-1] + chunkSize[ith];
	}
	thr.resize(num);
}

void Threads::
join()
{
	for (auto& o : thr) 
		o.join();
}

// Various chunksizes for flexible thread loads

inline void Threads::
set_chunks_equal( const szt w, const szt rest)
{
	for (szt ith=0; ith<num; ith++) {
		chunkSize[ith] = w/num;
		if (ith < rest) 
			chunkSize[ith]++;
	}
}

void Threads::
set_chunks_circular( const szt w, const szt rest)
{
// The coefficients are lengths of circle sagitta h = rnd * (1 - cos(phi/2)),
// where phi is the central angle (in radians)  defining the circle segment.
// The condition is the equality of the segment areas A = rnd^2 * ( phi - sin(phi) ) / 2, which reflect the thread loads.
// The problem has no explicit solution but can be solved numerically. E.g. in matlab:
//		% Let be given: 'k': an index in chunkSize[k] 
//		%				'num': the number of threads
//		% Then, to find 'h':
//		syms h;  vpasolve( 2*pi*k/num == 2*acos(1-h) - sin(2*acos(1-h)), h ) / 2
// 
	switch (num) {
	case 1 :											
		chunkSize[0] = w;	
		break;									 
	case 2 :											
		chunkSize[0] = szt(0.5000000 * w);						
		chunkSize[1] = w - chunkSize[0];						
		break;									 
	case 3 :											
		chunkSize[0] = szt(0.3675340 * w);
		chunkSize[1] = szt(0.6324660 * w) - chunkSize[0];
		chunkSize[2] =  w - chunkSize[1] - chunkSize[0];
		break;									 
	case 4 :											
		chunkSize[0] = szt(0.2980136 * w);
		chunkSize[1] = szt(0.5000000 * w) - chunkSize[0];
		chunkSize[2] = szt(0.7019864 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] =   w - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 5 :											
		chunkSize[0] = szt(0.2540691 * w);
		chunkSize[1] = szt(0.4211319 * w) - chunkSize[0];
		chunkSize[2] = szt(0.5788681 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.7459309 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] =   w - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 6 :											
		chunkSize[0] = szt(0.2233536 * w);
		chunkSize[1] = szt(0.3675340 * w) - chunkSize[0];
		chunkSize[2] = szt(0.5000000 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.6324660 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] = szt(0.7766464 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5] =   w - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 7 :											
		chunkSize[0] = szt(0.2004697 * w);
		chunkSize[1] = szt(0.3282611 * w) - chunkSize[0];
		chunkSize[2] = szt(0.4437815 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.5562185 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] = szt(0.6717389 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5] = szt(0.7995303 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6] =   w - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 8 :											
		chunkSize[0] = szt(0.1826477 * w);
		chunkSize[1] = szt(0.2980136 * w) - chunkSize[0];
		chunkSize[2] = szt(0.4011780 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.5000000 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] = szt(0.5988220 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5] = szt(0.7019864 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6] = szt(0.8173523 * w) - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[7] =   w - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 9 :											
		chunkSize[0] = szt(0.16830735 * w);
		chunkSize[1] = szt(0.27386929 * w) - chunkSize[0];
		chunkSize[2] = szt(0.36753396 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.45631111 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] = szt(0.54368889 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5] = szt(0.63246604 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6] = szt(0.72613071 * w) - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[7] = szt(0.83169265 * w) - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[8] =    w - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 10 :											
		chunkSize[0] = szt(0.15647559 * w);
		chunkSize[1] = szt(0.25406908 * w) - chunkSize[0];
		chunkSize[2] = szt(0.34015425 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3] = szt(0.42113190 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4] = szt(0.50000000 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5] = szt(0.57886810 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6] = szt(0.65984575 * w) - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[7] = szt(0.74593092 * w) - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[8] = szt(0.84352441 * w) - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[9] =    w - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 11 :											
		chunkSize[0]  = szt(0.14651773 * w);
		chunkSize[1]  = szt(0.23748417 * w) - chunkSize[0];
		chunkSize[2]  = szt(0.31735290 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3]  = szt(0.39205578 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4]  = szt(0.46426965 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5]  = szt(0.53573035 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6]  = szt(0.60794422 * w) - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[7]  = szt(0.68264710 * w) - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[8]  = szt(0.76251583 * w) - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[9]  = szt(0.85348227 * w) - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[10] =    w - chunkSize[9] - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	case 12 :											
		chunkSize[0]  = szt(0.13800066 * w);
		chunkSize[1]  = szt(0.22335364 * w) - chunkSize[0];
		chunkSize[2]  = szt(0.29801362 * w) - chunkSize[1] - chunkSize[0];
		chunkSize[3]  = szt(0.36753396 * w) - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[4]  = szt(0.43436113 * w) - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[5]  = szt(0.50000000 * w) - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[6]  = szt(0.56563887 * w) - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[7]  = szt(0.63246604 * w) - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[8]  = szt(0.70198638 * w) - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[9]  = szt(0.77664636 * w) - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[10] = szt(0.86199934 * w) - chunkSize[9] - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		chunkSize[11] =   w - chunkSize[10] - chunkSize[9] - chunkSize[8] - chunkSize[7] - chunkSize[6] - chunkSize[5] - chunkSize[4] - chunkSize[3] - chunkSize[2] - chunkSize[1] - chunkSize[0];
		break;									 
	default : 
		std::cout << "Error in Threads: Weights::CircleCenter is implemented for number of threads < 13 only. The attempted number is " << num << std::endl; 
		exit(0); 
	}
	
	for (szt ith=0; ith<rest; ith++) 
		chunkSize[ith]++;
}

void Threads::
set_chunks_triangleDecr( const szt size )
{
	szt w {};
	for (szt i=size; i>0; i--) 
		w += i;
	const szt e = w/num;
	const szt rest = w%num;

	std::vector<szt> dchs(num, e);
	for (szt i=0; i<rest; i++) 
		dchs[i]++;

	szt n(size);
	for (szt ith=0; ith<num-1; ith++) {
		chunkSize[ith] = 0;
		szt chs {};
		do {
			chunkSize[ith]++; 
			chs += n--;
		} while (chs < dchs[ith]);
	}
	chunkSize.back() = n;
}

void Threads::
print_regions( const bool withCout,
			   Msgr& msgr )
{
	if (withCout && msgr.so) {
		*msgr.so << " Thread borders: ";
		for (szt ith=0; ith<num; ith++) 
			*msgr.so << i1[ith] << " to " << i2[ith]-1 << "; ";
		*msgr.so << std::endl;
	}
	if (msgr.sl) {
		*msgr.sl << " Thread borders: ";
		for (szt ith=0; ith<num; ith++)
			*msgr.sl << i1[ith] << " to " << i2[ith]-1 << "; ";
		*msgr.sl << std::endl;
	}
}

}	// namespace Common
}	// namespace Utils
