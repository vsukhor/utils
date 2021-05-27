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

#include <fstream>
#include <iostream>
#include <random>
#include <type_traits>
#include <vector>

#include "../common/exceptions.h"
#include "../common/misc.h"
#include "../common/msgr.h"

/// Pseugo-random number generation.
namespace utils::random {

using szt = common::szt;
using Msgr = common::Msgr;
using szt = common::szt;
using uint = common::uint;

/// \brief Base class for random number factories.
/// \tparam realT Floating point type.
template <typename realT>
class Core {

    // Ensure that the template parameter is a floating type.
    static_assert(
        std::is_floating_point<realT>::value,
        "Class Core can only be instantiated with floating point types"
    );

public:

    ///< Number of seeds in the 'seed' file.
    static constexpr szt num_saved_seeds {1'000'001};
    ///< Size of the buffer for storing random numbers.
    static constexpr int bufferSize {1'000'000};
    ///< Master seed.
    static constexpr int mainSeed {1'234'567'890};

    /// \brief Produce \p num_saved_seeds seeds and store them in a file.
    /// \details This is done whenever the working directory does not already
    /// have such a file.
    static void make_seed(
        const std::string& seedFname,  ///< Name of the file with seeds.
        Msgr* msgr                     ///< Output message processor.
    );

    /// \brief Produce \p num_saved_seeds seeds and store them in a file.
    /// \details This is done whenever the working directory does not already
    /// have such a file.
    static uint readin_seed(
        const std::string& seedFname,  ///< Name of the file with seeds.
        szt runInd,                    ///< Run index to choose a seed.
        Msgr& msgr                     ///< Output message processor.
    );

    /// \brief \a seed getter.
    auto theSeed() { return seed; }

protected:

    explicit Core(
        Msgr& msgr,                   ///< Output message processor.
        uint seed,                    ///< Seed to use.
        const std::string& runName    ///< Human-readable run index.
    ) noexcept;

    explicit Core(
        Msgr& msgr,                    ///< Output message processor.
        const std::string& seedFname,  ///< Name of the file with seeds.
        szt runInd                     ///< Run index to choose a seed.
    );

private:

    uint  seed {common::huge<uint>};  ///< The seed
    Msgr& msgr;                       ///< Output message processor.
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename realT> 
Core<realT>::
Core(    Msgr& msgr,
        const uint seed,
        const std::string& runName) noexcept
    : seed {seed}
    , msgr {msgr}
{
    msgr.print("RUN " + runName);
    msgr.print("SEED = " + std::to_string(seed));
}

template <typename realT> 
Core<realT>::
Core(   Msgr& msgr,
        const std::string& seedFname,
        const szt runInd)
    : msgr {msgr}
{
    if (!common::file_exists(seedFname))
        make_seed(seedFname, &msgr);
    seed = readin_seed(seedFname, runInd, msgr);
    msgr.print("RUN = " + std::to_string(runInd));
    msgr.print("SEED = " + std::to_string(seed));
}

template <typename realT> 
void Core<realT>::
make_seed(
    const std::string& seedFname,
    Msgr* msgr
)
{
    if (const auto msg = "No seed file found. Creating a new seed file " + seedFname;
        msgr)
        msgr->print(msg);
    else
        std::cout << msg << std::endl;

    std::mt19937 g {mainSeed};

    constexpr int sd1 = 100'000'000;
    constexpr int sd2 = 2100'000'000;
    std::uniform_int_distribution<uint> seed_d(sd1, sd2);
    std::ofstream file {seedFname, std::ios::binary};
    if (!file.is_open()) {
        if (const auto msg = "Unable to create seed file "+seedFname;
            msgr)
            msgr->exit(msg);
        else
            common::exceptions::simple(msg, nullptr);
    }
    for (szt i=0; i<num_saved_seeds; i++) {
        const uint s = seed_d(g);
        file.write(reinterpret_cast<const char*>(&s), sizeof(uint));
    }
}

template <typename realT> 
uint Core<realT>::
readin_seed(
    const std::string& seedFname,
    szt runInd,
    Msgr& msgr
)
{
    msgr.print("Reading from file " + seedFname + " seed no: " +
               std::to_string(runInd));
    
    std::ifstream file {seedFname, std::ios::binary};
    if (!file.is_open())
        msgr.exit("Unable to open file "+seedFname);
    file.seekg(static_cast<std::fstream::off_type>(runInd * sizeof(uint)), file.beg);
    auto seed {common::huge<uint>};
    file.read(reinterpret_cast<char*>(&seed), sizeof(uint));

    return seed;
}

}    // namespace utils::random

#endif // UTILS_RANDOM_CORE_H
