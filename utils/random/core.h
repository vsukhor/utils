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

/// \file core.h
/// \brief Contains class Random::Core.
/// \author Valerii Sukhorukov

#ifndef UTILS_RANDOM_CORE_H
#define UTILS_RANDOM_CORE_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

#include "../common/misc.h"
#include "../msgr.h"

/// Pseugo-random number generation.
namespace utils::random {

using szt = common::szt;
using uint = common::uint;

/// \brief Base class for random number factories.
/// \tparam realT Floating point type.
template <typename realT>
class Core {

    // Ensure that the template parameter is a floating type.
    static_assert(
        std::is_floating_point_v<realT>,
        "Class Core can only be instantiated with floating point types"
    );

public:

    /// Size of the buffer for storing random numbers.
    static constexpr int bufferSize {1'000'000};
    /// Master seed.
    static constexpr int mainSeed {1'234'567'890};

    /// Produce \p num_saved_seeds seeds and store them in a file.
    /// This is done whenever the working directory does not already
    /// have such a file.
    static uint make_seed(
        uint runInd      ///< Run index to choose a seed.
    ) noexcept;

    /// Seed getter.
    auto get_seed() { return seed; }

protected:

    /// Constructor setting the seed uncoupled from run index.
    explicit Core(
        unsigned seed,                ///< Seed to use.
        const std::string& runName,   ///< Human-readable run index.
        Msgr& msgr                    ///< Output message processor.
    ) noexcept;

    /// Constructor setting the seed depending on run index.
    explicit Core(
        unsigned runInd,    ///< Run index to choose a seed.
        Msgr& msgr          ///< Output message processor.
    ) noexcept;

private:

    unsigned seed {common::huge<unsigned>};  ///< The seed
    Msgr& msgr;                              ///< Output message processor.
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename realT> 
Core<realT>::
Core(
    const unsigned seed,
    const std::string& runName,
    Msgr& msgr
) noexcept
    : seed {seed}
    , msgr {msgr}
{
    msgr.print("RUN = ", runName);
    msgr.print("SEED = ", seed);
}

template <typename realT> 
Core<realT>::
Core(
    const unsigned runInd,
    Msgr& msgr
) noexcept
    : seed {make_seed(runInd)}
    , msgr {msgr}
{
    msgr.print("RUN = ", runInd);
    msgr.print("SEED = ", seed);
}

template <typename realT> 
uint Core<realT>::
make_seed(
    const uint runInd
) noexcept
{
    std::mt19937 g {mainSeed};

    constexpr int sd1 = 100'000'000;
    constexpr int sd2 = 2'100'000'000;

    std::uniform_int_distribution<uint> seed_d(sd1, sd2);
    uint s {};
    for (uint i=0; i<runInd; i++)
        s = seed_d(g);

    return s;
}

}  // namespace utils::random

#endif // UTILS_RANDOM_CORE_H
