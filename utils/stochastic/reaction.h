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
 * \file reaction.h
 * \brief Interface for stochastically simulated reactions.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_STOCHASTIC_REACTION_H
#define UTILS_STOCHASTIC_REACTION_H

#include <string>
#include <utility>
#include <vector>

#include "../common/misc.h"
#include "../msgr.h"
#include "gillespie.h"

namespace utils::stochastic {

/// An abstract base class for all the reactions
template<typename RF>
class Reaction {

	friend Gillespie<RF, Reaction<RF>>;

public:  // Only constant parameters are public.

    using real = RF::Real;

	static std::vector<std::string>	activeInstances;

    /// Index in Simulation::rc, i.e. index among all used and not used reactions.
    const szt ind {undefined<szt>};

    const std::string shortName;  ///< Reaction name abbreviation.
    const std::string fullName;   ///< Extended reaction name.

	static bool is_active(
        const std::unique_ptr<Reaction>& r
    ) noexcept;

    /**
     * \brief Constructor.
     * \param [in,out] msgr Output message processor.
     * \param [in] ind Reaction id.
     * \param [in] rate Reaction rate constant.
     * \param [in] shortName Reaction name.
     * \param [in] fullName Reaction name abbreviated.
     */
    explicit Reaction(
        Msgr& msgr,
        szt ind,
        real rate,
        std::string shortName,  // value + move
        std::string fullName    // value + move
    );

	// The rule of five is triggered by the virtual destructor, the defaults suffice
    Reaction(const Reaction&) = default;             ///< Copy constructor.
    Reaction& operator=(const Reaction&) = default;  ///< Copy assignment.
    Reaction(Reaction&&) = default;                  ///< Move constructor.
    Reaction& operator=(Reaction&&) = default;       ///< Move assignment.
	virtual ~Reaction() noexcept = default; 		 ///< Destructor.

    /// Sets the Gillespie score for this reaction.
    virtual void set_score() noexcept = 0;
    /// Returns the Gillespie score for this reaction.
    real get_score() const noexcept { return *score; }

    /**
     * \brief Update propensity for a pair of network components.
     * \param [in] c1 Index of the 1st component to update.
     * \param [in] c2 Index of the 2nd component to update.
     */
    virtual void update_prop(szt c1, szt c2) noexcept = 0;
	virtual void update_prop(szt /* c1 */) noexcept {}
	virtual void set_prop() noexcept {}

    virtual void resize_pr() noexcept {}

    /// Execute the raction event.
    virtual void fire() noexcept = 0;

    /**
     * \brief Populates the vector of reactions that need a score update.
     * \details The update is performed after *this has fired
     *  Initializes the propensities and effective rate.
     * \param [in] rc Vector of unique pointers to all reactions taking part in
     *  the simulation.
     */
    virtual void initialize_dependencies(
        const std::vector<std::unique_ptr<Reaction>>& rc
    ) noexcept = 0;

    /**
     * \brief The number of times this reaction was fired.
     * \details \a eventCount getter.
     * \result The number of times this reaction was fired.
     */
    unsigned long event_count() const noexcept { return eventCount; }

    /**
     * \brief Prints the parameters common to all reactions.
     * \param le True if new line after the output.
     */
    virtual void print(bool le) const;

protected:

    Msgr& msgr;  ///< ref: Output message processor.

    /// Number of reaction events fired.
    unsigned long eventCount {};

    /// Reaction rate constant.
    real rate {};

    real* score {};

    /// Reactions that need a score update after *this has fired.
    std::vector<Reaction*> dependents;

    /**
     * \brief All necessary updates after the given reaction event was executed.
     * \details Pure virtual function: Network and reaction updates necessary
     *  after the given reaction event was executed.
     */
    virtual void update_netw_stats() = 0;

private:

    /**
     * \brief Attaches this score to the Gillespie mechanism.
     * \param [in] a Placeholder in the Gillespie object responsible for this
     *  reaction score.
     */
    void attach_score_pointer(real* a) noexcept { score = a; };
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename RF>
std::vector<std::string> Reaction<RF>::activeInstances;

template<typename RF>
inline
bool Reaction<RF>::
is_active(
    const std::unique_ptr<Reaction>& r
) noexcept
{
    return std::find(Reaction::activeInstances.begin(),
                     Reaction::activeInstances.end(),
                     r->shortName)
            != Reaction::activeInstances.end();
}

template<typename RF>
Reaction<RF>::
Reaction( Msgr& msgr,
          const szt ind,
          const real rate,
          const std::string shortName,  // value + move
          const std::string fullName    // value + move
    )
    : ind {ind}
    , shortName {std::move(shortName)}
    , fullName {std::move(fullName)}
    , msgr {msgr}
    , rate {rate}
{}

template<typename RF>
void Reaction<RF>::
print(const bool le) const
{
    msgr.print<false>(" shortName ", shortName);
    msgr.print<false>(" rate ", rate);
    msgr.print<false>(" score ", *score);
    msgr.print<false>(" eventCount ", eventCount);
    if (le) msgr.print("\n");
}

}  // namespace utils::stochastic

#endif  // UTILS_STOCHASTIC_REACTION_H
