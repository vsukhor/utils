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
* \file stop_watch.h
* Timing class.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_STOP_WATCH_H
#define UTILS_STOP_WATCH_H

#include <chrono>
#include <ctime>

#include "common/constants.h"

/// Library outer namespace.
namespace utils {

/// \struct StopWatch stop_watch.h
/// \brief Simple stop watch class using std::chrono::system_clock.
/// \details Implements convenience class for measuring time duration.
struct StopWatch {

    /// An instance in time.
    struct __attribute__((aligned(128)))
    Instance {

        /// A point in time.
        std::chrono::time_point<std::chrono::system_clock> h;

        std::time_t c {huge<std::time_t>};  ///< ctime-based instance.
        std::string str {};   ///< String-based representation.

        /// Current time.
        void operator()() {
            h = std::chrono::system_clock::now();
            c = std::chrono::system_clock::to_time_t(h);
            str = std::string(ctime(&c));
        }
    };

    Instance start;  ///< Start time.
    Instance stop;   ///< Stop time

    /// Duration between \a start and \a stop formatted as a string..
    std::string duration() {
        diff = stop.h - start.h;
        return STR(diff.count());
    }

private:

    /// Duration in std::chrono format.
    std::chrono::duration<double> diff {huge<double>};

};

}  // namespace utils

#endif // UTILS_STOP_WATCH_H
