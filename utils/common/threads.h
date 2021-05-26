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
* \file threads.h
* Contains class Threads.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_COMMON_THREADS_H
#define UTILS_COMMON_THREADS_H

#include <iostream>
#include <thread>
#include <vector>

#include "misc.h"
#include "msgr.h"

/// General stuff.
namespace Utils::Common {

/**
* \class Threads threads.h
* \brief Simple threading utility class.
* \details Implements convenience class for handling a collection
* of std::thread objects.
*/
class Threads {

public:

    /**
    * \brief Enumerates basic load sharing modes.
    * \details Names three modes of load distribution between threads.
    */
    enum class Weights {
        CircleCenter,   ///< Circle-shaped load distribution.
        Equal,          ///< Uniform load sharing.
        TriangleDecr    ///< Right triangle-shaped distribution
    };

    const szt        num;     ///< Number of threads.
    std::vector<szt> chs;     ///< Per-thread amounts of relative load.
    std::vector<szt> i1, i2;  ///< Range borders.

    std::vector<std::thread> thr;  ///< Container holding the threads.

    /**
    * \brief Constructor.
    * \details Creates threads based on a set of work units.
    * \param offset Offset from the start of work unit container.
    * \param size Size of the work unit container shared among the threads.
    * \param omittedBoundaries Boundsary work units to discard.
    * \param wht Relative weiting.
    * \param nThreads Thread number.
    */
    explicit Threads(
        szt offset,
        szt size,
        ulong omittedBoundaries,
        const Weights wht,
        ulong nThreads );
    
    /**
    * Joins the threads.
    */
    void join();
    
    // Various weights for relative thread loads
    /**
    * \brief Sets weighting factors according to \a Weights::Equal.
    * \param w Total number of work units.
    * \param rest Number of work units remaining after the optimal distribution.
    */
    void set_chunks_equal(szt w, szt rest);

    /**
    * \brief Sets weighting factors according to \a Weights::CircleCenter.
    * \param w Total number of work units.
    * \param rest Number of work units remaining after the optimal distribution.
    */
    void set_chunks_circular(szt w, szt rest);
    /**
    * \brief Sets weighting factors according to \a Weights::TriangleDecr.
    * \param size Total number of work units.
    */
    void set_chunks_triangleDecr(szt size);

    /**
    * \brief Prints work unit borders for particular threads.
    * \param withCout Specifies if printing to cout.
    * \param msgr \a Msgr used for the output.
    */
    void print_regions(bool withCout,
                       Msgr& msgr);
};

}  // namespace Utils::Common

#endif // UTILS_COMMON_THREADS_H
