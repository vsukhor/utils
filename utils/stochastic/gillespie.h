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
* \file gillespie.h
* \brief Contains class Gillespie.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_STOCHASTIC_GILLESPIE_H
#define UTILS_STOCHASTIC_GILLESPIE_H

#include <algorithm>
#include <cmath>
#include <fstream>

#include "../common/constants.h"
#include "../common/misc.h"

/// General stuff.
namespace utils::stochastic {

/**
 * \class Gillespie gillespie.h
 * \brief Exact Stochastic Simulation algorithm.
 * \details Implementation of the Exact Stochastic Simulation algorithm
 * of D.T. Gillespie <https://pubs.acs.org/doi/abs/10.1021/j100540a008>
 * \tparam RF random factory class.
 * \tparam Reaction base class for reactions.
*/
template <typename RF, typename Reaction>
class Gillespie {

public :

    using realT = float;
    using szt = common::szt;
    template <typename T> using vup = std::vector<std::unique_ptr<T>>;
    static constexpr auto hugeszt = common::huge<szt>;
    static constexpr auto real0 = static_cast<realT>(0);
    static constexpr auto real1 = static_cast<realT>(1);
    /**
     * \brief Constructor.
     * \param rnd Random class factory.
    */
    explicit Gillespie(RF& rnd) noexcept;

    /**
     * \brief Adds a reaction \p r to the reaction collection.
     * \param r Unique pointer to the new reaction.
    */
    void add_reaction(std::unique_ptr<Reaction> r);

    /**
     * \brief Performs neccessary settings of the class members.
    */
    void initialize() noexcept;

    bool set_asum() noexcept;
    void fire(double& time) noexcept;

    /**
     * \brief Getter for the time till the nest reaction event.
     * \return Time till the nest reaction event.
    */
    [[nodiscard]] realT tau() const noexcept { return tau_; }

    /**
     * \brief Prints reaction propensities to \p os.
     * \param os Stream to print to.
    */
    void printScores(std::ostream& os) const;

    /**
     * \brief Prints a log record with the information on current reaction.
     * \param os Stream to record to.
    */
    void log_data(std::ostream& os) const;

    /**
     * \brief Number of reactions.
     * \return Number of reactions.
    */
    [[nodiscard]] constexpr auto num_reactions() const noexcept -> szt;

    /**
     * \brief Short human readable name of the reaction.
     * \param ind index of the reaction.
     * \returns Short human readable name of the reaction indexed by \p ind.
    */
    [[nodiscard]] std::string short_name(szt ind) const noexcept;

    /**
     * \brief Finds reaction from its name.
     * \param name Name of the reaction to find.
     * \returns Pointer to the reaction \p name ot nullptr if does not exist.
    */
    Reaction* get_reaction(const std::string& name) const;

private:

    vup<Reaction>      rc;         ///< Vector of unique pointers to reactions.
    std::vector<realT> a;          ///< Vector of protensitioe.
    szt                rind {hugeszt}; ///< Index of the current reaction.
    realT              tau_ {};    ///< Time till the next reaction event.
    RF&                rnd;        ///< Random number generator.
    szt                nreact {};  ///< Total number of reactions.
    std::vector<szt>   rtype;      ///< Human-readable abbrev. of reaction names.
    std::vector<realT> csums;      ///< Cumulative sum.
    realT              asum {};    ///< Total propensity of the reaction set.
    std::vector<realT> auxf;       ///< Auxiliary vector of nreact reals.
    std::vector<szt>   auxi;       ///< Auxiliary vector of nreact integers.
    std::vector<szt>   rinds;      ///< Reaction indexes.

    /**
     * \brief Sets the time \tau_ till the next reaction event.
    */
    void set_tau() noexcept;

    /**
     * \brief Sets index of the reaction to fire in the next event.
    */
    void set_rind() noexcept;

    /**
     * \brief Produces pointer to the current reaction.
     * \return Pointer to the current reaction.
    */
    constexpr Reaction* currRc() const noexcept;

};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename RF, typename Reaction> inline
Gillespie<RF,Reaction>::
Gillespie( RF& rnd ) noexcept
    : rnd {rnd}
{}

template <typename RF, typename Reaction> inline
void Gillespie<RF,Reaction>::
add_reaction( std::unique_ptr<Reaction> rup )
{
    rc.push_back(std::move(rup));
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
initialize() noexcept
{
    nreact = rc.size();
    for (szt ri=0; ri<nreact; ri++)
        rtype.push_back(ri);

    a.resize(nreact);
    for (szt i=0; i<nreact; i++) {
        rc[rtype[i]]->attach_score_pointer(&a[i]);
        rc[rtype[i]]->initialize_dependencies(rc);
    }
    
    auxf.resize(nreact);
    csums.resize(nreact);
    auxi.resize(nreact);
    rinds.resize(nreact);
}

template <typename RF, typename Reaction>
Reaction* Gillespie<RF,Reaction>::
get_reaction( const std::string& name ) const
{
    for (const auto& o : rc)
        if (o->shortName == name)
            return o.get();

    return nullptr;
}

template <typename RF, typename Reaction>
bool Gillespie<RF,Reaction>::
set_asum() noexcept
{
    asum = std::accumulate(a.begin(), a.end(), real0);
    
    return asum != real0;
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
set_rind() noexcept
{
    for (szt i=0; i<nreact; i++)
        auxf[i] = a[i] / asum;
    std::partial_sum(auxf.begin(), auxf.end(), csums.begin());

    realT ran {};
    do ran = rnd.r01u(); 
    while (ran >= csums[nreact-1]);

    for (szt i=0; i<nreact; i++)
        auxi[i] = szt(ran < csums[i]);

    const szt rindnum {common::find(auxi, rinds)};
    rind = *std::min_element(rinds.begin(), rinds.begin()+rindnum);
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
set_tau() noexcept
{
    realT ran {};
    do ran = rnd.r01u(); 
    while (ran <= real0 || ran >= real1);

    tau_ = std::log(real1 / ran) / asum;
    XASSERT(!std::isnan(tau_), "Tau is nan");
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
fire( double& time ) noexcept
{
    set_rind();
    set_tau();
    time += tau_;
//    (*currRc())();
    currRc()->fire();
}

template <typename RF, typename Reaction> constexpr
Reaction* Gillespie<RF,Reaction>::
currRc() const noexcept
{
    return rind < hugeszt ? rc[rtype[rind]].get()
                            : nullptr;
}

template <typename RF, typename Reaction> constexpr
auto Gillespie<RF,Reaction>::
num_reactions() const noexcept -> szt
{
    return rc.size();
}

template <typename RF, typename Reaction>
std::string Gillespie<RF,Reaction>::
short_name( const szt i ) const noexcept
{
    return rc[i]->shortName;
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
log_data( std::ostream& os ) const
{
    auto pad_zeros = [](auto n) {
        constexpr int MIN2 = 10;
        constexpr int MIN3 = 100;
        return (n < MIN2 ? "00"
                         : (n < MIN3 ? "0"
                                     : "")) + std::to_string(n);
    };

    os << " tau " << tau_
       << " rt " << (rind!=hugeszt ? "" : "000")
                  << pad_zeros(rind);
    if (rind != hugeszt)
        os << " " << currRc()->shortName;
}

template <typename RF, typename Reaction>
void Gillespie<RF,Reaction>::
printScores( std::ostream& os ) const
{
    for (szt i=0; i<nreact; i++)
        rc[i]->set_score();
    os << " scores: ";
    for (szt i=0; i<nreact; i++)
        os << rc[i]->shortName << " "
           << rc[i]->event_count() << " "
           << rc[i]->get_score() << " ";
}

}  // namespace utils::stochastic

#endif // UTILS_STOCHASTIC_GILLESPIE_H
