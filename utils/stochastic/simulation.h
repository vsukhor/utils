/* =============================================================================

 Copyright (C) 2009-2025 Valerii Sukhorukov. All Rights Reserved.

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
 * \file simulation.h
 * \brief Common-ground interface for simulation classes.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_STOCHASTIC_SIMULATION_H
#define UTILS_STOCHASTIC_SIMULATION_H

#include <cstdint>
#include <ostream>
#include <string>

#include "../msgr.h"
#include "../constants.h"
#include "gillespie.h"
#include "reaction.h"

namespace utils::stochastic {

using itT = std::uint_fast64_t; ///< Type for counting simulation iterations.

/**
 * \brief Interface defining essential components of simulation classes.
 * \tparam System Typename of the simulated system.
 * \tparam Engine Typename of the engine controlling the simulation progrression.
 */
template<typename System,
         typename Engine>
class Simulation {

public:

    explicit Simulation(
        System& obj,
        const std::string& runName,
        Engine& engine,
        double& time,
        itT& it,
        Msgr& msgr
    );

    virtual ~Simulation() = default;

protected:

    using ReactionBase = Reaction<typename Engine::RandFactory>;

    System& obj;  ///< Ref to the simulated system.

    Msgr& msgr;  ///< Logging.

    /// Human-readable name of the simulation run.
    const std::string& runName;

    /// Stochastic engine controlling the simulation progrression.
    Engine& engine;

    double& time;  ///< Ref to simulation time.
    itT&    it;    ///< Ref to iteration count

    itT logcnt {};  ///< Counter for log reports.

    /// Flag marking if the system is currently considered equilibrated.
    bool isEquilibrated {};

    // The rule of five is triggered by the destructor, the defaults suffice:
    Simulation(const Simulation&) = default;             // Copy constructor.
    Simulation& operator=(const Simulation&) = default;  // Copy assignment.
    Simulation(Simulation&&) = default;                  // Move constructor.
    Simulation& operator=(Simulation&&) = default;       // Move assignment.

    /// Main loop running the simulation in continuous time.
    virtual void run_tcont() = 0;

    /// Handler for the case when all reactions are exosted.
    void asum0fun();

    /// Output a log line reporting current system state to a specific stream.
    virtual void log(std::ostream&) = 0;

    /// Output a log line reporting current system state to all Msgr streams.
    void log();
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename System,
         typename Engine>
Simulation<System, Engine>::
Simulation(
    System& obj,
    const std::string& runName,
    Engine& engine,
    double& time,
    itT& it,
    Msgr& msgr
)
    : obj {obj}
    , msgr {msgr}
    , runName {runName}
    , engine {engine}
    , time {time}
    , it {it}
{}

template<typename System,
         typename Engine>
void Simulation<System, Engine>::
asum0fun()
{
    msgr.print("\n");
    log();
    msgr.print("No reaction left! ",
               "Termination due to reaction *score == 0 for all reactions used.");
}

template<typename System,
         typename Engine>
void Simulation<System, Engine>::
log()
{
    if (msgr.so) log(*msgr.so);
    if (msgr.sl) log(*msgr.sl);
    logcnt++;
}

}  // namespace utils::stochastic

#endif  // UTILS_STOCHASTIC_SIMULATION_H
