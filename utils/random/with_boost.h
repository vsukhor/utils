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
* \file with_boost.h
* Contains class Boost.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_RANDOM_WITH_BOOST_H
#define UTILS_RANDOM_WITH_BOOST_H

#include <cmath>
#include <filesystem>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/random/binomial_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "../arrays/all.h"
#include "../common/constants.h"
#include "../common/misc.h"
#include "core.h"

/// Pseugo-random number generation.
namespace utils::random {

/// \brief Random number factory based on Boost distribution functions.
/// \tparam realT Floating point type.
template <typename realT>
class Boost
    : public Core<realT> {

    // Ensure that the template parameter is a floating type
    static_assert(
        std::is_floating_point<realT>::value,
        "Class Core can only be instantiated with floating point types"
    );


public:

    using real = realT;
    using A2r = arrays::A2<realT>;
    using A3r = arrays::A3<realT>;
    static constexpr auto half = common::half<realT>;
    static constexpr auto one = common::one<realT>;
    static constexpr auto two = common::two<realT>;
    static constexpr auto zero = common::zero<realT>;
    static constexpr auto pi = common::pi<realT>;
    static constexpr auto twopi = common::twopi<realT>;
    static constexpr auto halfpi = common::halfpi<realT>;

//    /// \brief Default constructor.
//    Boost() = default;

    /// \brief Constructor.
    /// \param seedFile Name of the file contining seeds.
    /// \param ii Run index.
    /// \param msgr Output message processor.
    explicit Boost(
        const std::filesystem::path& seedFile,
        szt ii,
        Msgr& msgr);

    /// \brief Constructor.
    /// \param seed Random number generator seed.
    /// \param runName Human-readable run index.
    /// \param msgr Output message processor.
    explicit Boost(
        uint seed,
        const std::string& runName,
        Msgr& msgr
    );

    /// A pseudo-random number with uniform distribution over [0,1).
    realT r01u();
    
    /// \brief A pseudo-random signed int from the range [0, max-1].
    /// \param max Max boundary of the sampled range.
    constexpr int uniform0(int max);

    /// \brief A pseudo-random unsigned int from the range [0, max-1].
    /// \param max Max boundary of the sampled range.
    constexpr uint uniform0(uint max);

    /// \brief A pseudo-random std::size_t from the range [0, max-1].
    /// \param max Max boundary of the sampled range.
    constexpr szt uniform0(szt max);

    // A pseudo-random realT from the range [0., max].
    /// \param max Max boundary of the sampled range.
    constexpr realT uniform0(realT max);
    
    /// \brief A pseudo-random integer from the range [1, max].
    /// \tparam intT Integer fundamental type.
    /// \param max Max boundary of the sampled range.
    template <typename intT>
    constexpr intT uniform1(intT max);
    
    /// \brief A point uinformly distributed within \p solidAngle on the surface of a unit sphere.
    /// \param solidAngle Constraining solid angle.
    constexpr auto uniform_direction(
        realT solidAngle  //=pi
    ) noexcept -> A3r;

    /// \brief A point uinformly distributed within boundaries on the surface of a unit sphere.
    /// Inclination is limited by \p inclMinMax [0, pi) around +z axis
    /// direction (i.e. \p phPole == 0).
    /// Azimuth is limited by \p azimMinMax [-pi, pi) around +x axis
    /// direction  (i.e. \p th == 0).
    /// \param[in] inclMinMax Min, max constrains on inclination.
    /// \param[in] azimMinMax Min, max constrains on azimuth.
    /// \param[in] azimSymmetric Switch if the szimuth is symmetric.
    /// \param[out] phPole Inclination of the resulting point.
    /// \param[out] th Azimuth of the resulting point.
    /// \return Point on a sphere uinformly distributed within angular boundaries.
    constexpr auto uniform_direction(
        const A2r& inclMinMax,
        const A2r& azimMinMax,
        bool azimSymmetric,
        realT& phPole,
        realT& th
    ) noexcept -> A3r;

    /// A point uinformly distributed within boundaries on a sphere surface.
    /// Implements trigonometric method.
    /// \param solidAngle Surface patch where the random point may belong to;
    /// set it to pi for the whole surface.
    /// \param r Shpere of radius.
    /// \param poleDir [0,1,2] is the index of the axis around which \p
    /// solidAngle is set.
    /// \return A point uinformly distributed within \p solidAngle on a
    /// shpere of radius \p r.
    auto uniform_on_sphere(
        realT solidAngle,  //=pi,
        realT r,  //=one,
        int poleDir    //=2
    ) noexcept -> A3r;
    
    /// \brief A point uinformly distributed within boundaries on a spheroid surface.
    /// Implements trigonometric method.
    /// \param solidAngle Surface patch where the random point may belong to;
    /// set it to pi for the whole surface.
    /// \param r Spheroid semi-axes dimensions: r[0] = a = b, r[1] = c
    /// \param poleDir [0,1,2] is the index of the axis around which \p
    /// solidAngle is set.
    /// \param bias [-1,0,1]: -1 towards poles; 1 towards equator, 0 no bias
    /// \param biasPar Bias strength.
    /// \return A point uinformly distributed within \p solidAngle on a
    /// shperoid of semi-axes dimensions \p r.
    auto uniform_on_spheriod(
        realT solidAngle, //=pi,
        const A2r& r,  //=one,
        int poleDir,  //=2,
        int bias,  //=0,
        realT biasPar  //=zero
    ) noexcept -> A3r;

    /// \brief Uinform directional (angular) distribution
    /// \details A point on the ellipse boundary from a uinform directional
    /// (angular) distribution.
    /// \note This is not a uniform density over the ellipse boundary.
    /// The ellipse is centered at (0,0)
    /// \param r Semi-axes dimensions of the ellipse: r = {a, b}.
    constexpr auto uniform_on_ellipse(
        const A2r& r  //=one
    ) noexcept -> A2r;

    /// \brief Point uinformly distributed over ellipse area .
    /// \details A point within the ellipse boundary having uinform
    /// distribution over the ellipse area.
    /// The ellipse is centered at (0,0)
    /// \param r Semi-axes dimensions of the ellipse: r = {a, b}.
    constexpr auto uniform_in_ellipse(
        const A2r& r  //=one
    ) noexcept -> A2r;

/*    /// \brief A shifted point within the ellipse boundary having uinform
    /// distribution over the ellipse area.
    /// \param r1 Semi-axes dimensions of the ellipse: r1 = {a, b}.
    /// \param r0 Ellipse center.
    /// \param shift Shift.
    A2<realT> uniform_in_ellipse(const A2<realT>& r1,
                                 const A2<realT>& r0, 
                                 const A2<realT>& shift=zero) noexcept;
*/
    /// \brief Exponentially distributed random numbers.
    /// \return A pseudo-random number sampled from exponential distribution.
    constexpr realT exponential_number(
        realT mi      ///< Rate parameter.
    ) noexcept;

    /// \brief Poisson distributed random numbers.
    /// \return A pseudo-random number sampled from Poisson distribution.
    uint poisson_number(
        realT lambda  ///< Rate parameter.
    ) noexcept;

    /// \brief Weibull distributed random numbers.
    /// \return A pseudo-random number sampled from Weibull distribution.
    constexpr realT weibull_number(
        realT lambda,  ///< Scale parameter.
        realT k        ///< Shape parameter.
    )  noexcept;

    /// \brief Logistically distributed random numbers.
    /// \return A pseudo-random number sampled from logistic distribution.
    constexpr realT logistic_number(
        realT mi,      ///< Mean.
        realT s        ///< Scale parameter.
    )  noexcept;

    /// \brief Binomially distributed (n, p) pseudo-random number.
    /// \return A pseudo-random number sampled from binomial distribution.
    uint binomial_number(
        uint n,    ///< Number of trials.
        realT p    ///< Outcome probability.
    )  noexcept;

    /// \brief Multinomially distributed pseudo-random numbers.
    /// \details Of the \p n independent trials each of which leads to a
    /// success for exactly one of \p k categories,
    /// with each category having a given fixed success probability.
    /// \return A vector of multinomially distributed deviates.
    std::vector<uint> multinomial_number(
        uint n,   ///< Number of trials.
        uint k    ///< Number of categories.
    )  noexcept;

    /// \brief Multinomially distributed pseudo-random numbers.
    /// \details Distributes \p n into p.size()+1 slots with
    /// probabilities p[0], p[1], ..., p.back(), 1 - sum(p)
    /// \return A vector of multinomially distributed deviates.
    std::vector<uint> multinomial_number(
        uint n,               ///< Number of trials.
        std::vector<realT> p  ///< Vector of probabilities.
    ) noexcept;

    /// \brief Normally distributed pseudo-random number.
    /// \return Normally distributed deviate N(mi, sigma^2).
    constexpr realT gaussian_number(
        realT mi,     ///< Mean.
        realT sigma   ///< Standard deviation.
    ) noexcept;

    /// \brief Constrained normal deviate.
    /// \details Normally distributed pseudo-random number constrained
    /// between \p cmin and \p cmax.
    /// \return Normally distributed deviate N(mi, sigma^2).
    constexpr realT gaussian_number_constrained(
        realT mi,     ///< Mean.
        realT sigma,  ///< Standard deviation.
        realT cmin,   ///< Min boundary.
        realT cmax    ///< Max boundary.
    ) noexcept;
    
private:

    using Core<realT>::bufferSize;

    /// Uniform 0 to 1 float.
    boost::random::uniform_01<float>  flt01_unifromDistr;
    /// Uniform 0 to 1 double.
    boost::random::uniform_01<double> dbl01_unifromDistr;
    /// Standard normal distribution.
    boost::random::normal_distribution<realT> normalDistr;

    /// Buffer array for storing random numbers.
    std::array<realT, bufferSize> rU01;

    /// Index of the current random number in \a rU01.
    volatile long rU01_ind;

//    boost::mt19937 g;   ///< Random number generator.
    std::mt19937 g;       ///< Random number generator.

    /// Populate the buffer array \a rU01 with a new butch of random numbers.
    void prepare_uniform_real01();
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename realT>
Boost<realT>::
Boost( const std::filesystem::path& seedFile,
       const szt ii,
       Msgr& msgr)
    : Core<realT> {msgr, seedFile, ii}
    , rU01_ind {-1}
{
    g.seed(this->theSeed());
    prepare_uniform_real01();
}


template <typename realT>
Boost<realT>::
Boost( const uint seed,
       const std::string& runName,
       Msgr& msgr)
    : Core<realT> {msgr, seed, runName}
    , rU01_ind {-1}
{
    g.seed(this->theSeed());
    prepare_uniform_real01();
}


// Generates realT random numbers with uniform distribution over [0,1)
template <> inline
void Boost<float>::
prepare_uniform_real01()
{
    for (auto& o : rU01) 
        o = flt01_unifromDistr(g);    
}


template <> inline
void Boost<double>::
prepare_uniform_real01()
{
    for (auto& o : rU01) 
        o = dbl01_unifromDistr(g);    
}


// Returns a random number with uniform distribution over [0,1).
template <typename realT> inline
realT Boost<realT>::
r01u()
{
    auto counter = rU01_ind + 1;  // local because of rU01_ind volatility
    if (counter == bufferSize) {
        prepare_uniform_real01();
        rU01_ind = 0;
    }
    else rU01_ind = counter;

    return rU01[rU01_ind];
}


// Returns int in the range [0, max-1]
template <typename realT> constexpr
int Boost<realT>::
uniform0( const int max )
{            
    XASSERT(max > 0, "Boost<realT>::uniform0 requires max > 0");

    auto ir {static_cast<int>(r01u() * max)};
    
    while (ir >= max) 
        ir = static_cast<int>(r01u() * max);
    
    return ir;    
}


// Returns uint in the range [0, max-1].
template <typename realT> constexpr
uint Boost<realT>::
uniform0( const uint max )
{            
    XASSERT(max > 0, "Boost<realT>::uniform0 requires max > 0");

    auto ir = static_cast<uint>(r01u() * max);
    
    while (ir >= max) 
        ir = static_cast<uint>(r01u() * max);
    
    return ir;    
}


// Returns szt in the range [0, max-1].
template <typename realT> constexpr
szt Boost<realT>::
uniform0(const szt max)
{            
    XASSERT(max > 0, "Boost<realT>::uniform0 requires max > 0");

    auto ir = static_cast<szt>(r01u() * max);
    
    while (ir >= max) 
        ir = static_cast<szt>(r01u() * max);
    
    return ir;    
}


// Returns outT in the range [1, max].
template <typename realT>
template <typename intT> constexpr
intT Boost<realT>::
uniform1( const intT max )
{            
    // Ensure that the template parameter is a floating type
    static_assert(std::is_integral<intT>::value,
                  "This function can only be instantiated with integer types");
    
    XASSERT(max > 0, "Boost<realT>::uniform1 requires max > 0 ");
    
    return uniform0(max) + 1;
}

template <typename realT> constexpr
realT Boost<realT>::
uniform0(const realT max)
{            
    XASSERT(max > zero, "Boost<realT>::uniform0 requires max > 0 ");

    auto ir = r01u() * max;    
    
    while (ir >= max) 
        ir = r01u() * max;
    
    return ir;    
}


// Returns a point uinformly distributed within solidAngle on a shpere.
template <typename realT> constexpr
auto Boost<realT>::
uniform_direction( const realT solidAngle ) noexcept -> A3r
{
    do {
        // inclination of the candidate point on a sphere surface
        const auto ph = pi * (r01u() - half);
        if (ph > solidAngle) continue;
        
        // azimuth of the candidate point on a sphere surface
        auto th = twopi * (r01u() - half);
        
        // if upper altitude paralleles are shorter, reject some
        // points positioned outside their length:
        if (std::abs(th) < pi * std::cos(ph)) {
            // Spread the remaining points over the length of the parallele:
            th /= std::cos(ph);
            
            return { std::cos(ph) * std::cos(th),
                     std::cos(ph) * std::sin(th),
                     std::sin(ph) };
        }
    } while (true);
    
    return {};
}


// Returns a point uinformly distributed on a shpere within inclMinMax
// and azimMinMax. Inclination is limited by inclMinMax [0, pi) around +z
// axis direction (i.e. phPole == 0). Azimuth is limited by azimMinMax [-pi, pi)
//  around +x axis direction  (i.e. th == 0).
template <typename realT> constexpr
auto Boost<realT>::
uniform_direction(
    const A2r& inclMinMax,
    const A2r& azimMinMax,
    const bool azimSymmetric,
    realT& phPole,
    realT& th
) noexcept -> A3r
{
    while (true) {
        // Inclination of the candidate point on a sphere surface:
        const auto ph = pi*r01u() - halfpi;
        if (ph >  halfpi - inclMinMax[0] ||
            ph <= halfpi - inclMinMax[1]) continue;

        // Azimuth of the candidate point on a sphere surface:
        th = twopi * r01u() - pi;    // [-pi, pi)

        // If upper altitude paralleles are shorter, reject some points
        // positioned outside their length:
        if (std::abs(th) < pi * std::cos(ph)) {
            // 'th' is spread the remaining points over the length of the parallele.
            th /= std::cos(ph);
            bool reject = th <  azimMinMax[0] ||
                          th >= azimMinMax[1];
            if (azimSymmetric)
                reject = reject &&
                        (th >= -pi + azimMinMax[1] &&
                         th <   pi + azimMinMax[0]);
            if (reject) continue;

            phPole = halfpi - ph;    // relative to the pole

            return { std::cos(ph) * std::cos(th),
                     std::cos(ph) * std::sin(th),
                     std::sin(ph) };
        }
    }

    return zero;
}


// Trigonometric method: returns a point uinformly distributed within
// solidAngle on a shpere of radius r.
// 'solidAngle' surface patch where the random point may belong to; set it
//      to pi for the whole surface.
// 'poleDir' [0,1,2] is the direction of the solidAngle axis.
template <typename realT> 
auto Boost<realT>::
uniform_on_sphere(
    const realT solidAngle,
    const realT r,
    const int poleDir
) noexcept -> A3r
{    
    XASSERT(solidAngle > zero &&
            solidAngle <= pi,
            "Error in Random::uniform_on_sphere: incorrect solidAngle");
    XASSERT(r > zero,
            "Error in Random::uniform_on_sphere: incorrect r");
    XASSERT(poleDir >= 0 &&
            poleDir <= 2,
            "Error in Random::uniform_on_sphere: incorrect poleDir");

    while (true) {
        const auto u = two * r01u() - one;
        if (u < std::cos(solidAngle))
            // reject the points outside the cap set by solidAngle,
            // i.e one with hieght h = r * (1 - std::cos(solidAngle))
            continue;
        
        const auto v = std::sqrt(one - u*u);
        // Inclination of the candidate point on a sphere surface:
        const auto ph = twopi * r01u();

        const auto vcp = v * std::cos(ph);
        const auto vsp = v * std::sin(ph);

        return poleDir == 2
               ? A3r(vcp, vsp, u) * r
               : (poleDir == 0
                  ? A3r(u,   vcp, vsp) * r
                  : A3r(vcp, u,   vsp) * r);
    }
}


// Trigonometric method: returns a point uinformly distributed within solidAngle
// on a shperoid of dimensions r;
// 'solidAngle' surface patch where the random point may belong to;
// set it to pi for the whole surface.
// 'r' is spheroid dimensions: r[0] = a = b, r[1] = c
// 'poleDir' [0,1,2] is the direction of the solidAngle axis
// 'bias' [-1,0,1]: -1 to poles; 1 to equator, 0 none
template <typename realT> 
auto Boost<realT>::
uniform_on_spheriod(
    const realT solidAngle,
    const A2r& r,
    const int poleDir,
    const int bias,
    const realT biasPar
) noexcept -> A3r
{
    // Check that the pole direction is along one of the major axes: 0, 1, 2:
    XASSERT(solidAngle > zero &&
            solidAngle <= pi,
            "Error in Random::uniform_on_spheriod: incorrect solidAngle");
    XASSERT(r > zero,
            "Error in Random::uniform_on_spheriod: incorrect r");
    XASSERT(poleDir >= 0 &&
            poleDir <= 2,
            "Error in Random::uniform_on_spheriod: incorrect poleDir");
    XASSERT(!bias || bias == -1 || bias == 1,
            "Error in Random::uniform_on_spheriod: incorrect bias");

    while (true) {
        // Inclination of the candidate point on a sphere surface:
        const auto phi = twopi * r01u();

        const auto u = two * r01u() - one;
        const auto v = std::sqrt(one - u*u);

        const auto x = r[0] * v * std::cos(phi);
        const auto y = r[0] * v * std::sin(phi);
        const auto z = r[1] * u;

        // Reject the points outside the cap set by solidAngle,
        // i.e one with hieght h = r * (1 - std::cos(solidAngle)):
        if (const auto sathr = std::cos(solidAngle);
            (poleDir == 2 && z < r[1] * sathr) ||
            (poleDir == 1 && y < r[0] * sathr) ||
            (poleDir == 0 && x < r[0] * sathr))
            continue;
        
        const auto s0 = std::sqrt(
            (x*x + y*y) / r[0]*r[0]*r[0]*r[0] +
                    z*z / r[1]*r[1]*r[1]*r[1]
        );
         // Acceptance threshold for prolate vs. oblate case:
        if (const auto s = s0 * (r[0] < r[1] ? r[0] : r[1]);
            s < r01u())
            continue;

        const A3r res {x, y, z};

        if (bias == -1)  // Bias towards poles.
            if (std::acos(std::abs(z) / res.norm()) >
                gaussian_number_constrained(zero, biasPar, zero, halfpi))
                continue;

        if (bias == 1) // Bias towards equator.
            if (std::acos(std::sqrt(x*x + y*y) / res.norm()) >
                 gaussian_number_constrained(zero, biasPar, zero, halfpi))
                continue;
        // No bias.
        return res;
    } 
}


template <typename realT> constexpr
auto Boost<realT>::
uniform_on_ellipse(
    const A2r& r    ///< ellipse dimensions: r = {a, b}
) noexcept -> A2r
{
    // Inclination of the candidate point:
    const auto phi = twopi * r01u();

    return { r[0] * std::cos(phi),
             r[1] * std::sin(phi) };

}


template <typename realT> constexpr
auto Boost<realT>::
uniform_in_ellipse(
    const A2r& r    ///< ellipse dimensions: r = {a, b}
) noexcept -> A2r
{
    const auto rho = r01u();
    // Inclination of the candidate point:
    const auto phi = twopi * r01u();
    
    // (x, y) is a random point inside a circle of radius 1
    const auto x = std::sqrt(rho) * std::cos(phi);
    const auto y = std::sqrt(rho) * std::sin(phi);
    
    return {x * r[0],
            y * r[1]};
}
/*
template <typename realT> constexpr
A2<realT> Boost<realT>::
uniform_in_ellipse( const A2<realT>& r1,
                    const A2<realT>& r0,
                    const A2<realT>& shift ) noexcept
{    
    do
        if (const auto p = uniform_in_ellipse(r1);
            p[0] >= shift[0] + r0[0] ||
            p[1] >= shift[1] + r0[1]) 
            return p;
    while (true);
}
*/


template <typename realT> constexpr
realT Boost<realT>::
exponential_number(
    const realT mi
) noexcept
{
    XASSERT(mi >= zero, "Boost<realT>::exponentialNum requires mi >= 0");

    return - mi * std::log(r01u());
}


// Returns a poissonian distributed deviate with mean mi.
template <typename realT>
uint Boost<realT>::
poisson_number(
    const realT lambda
) noexcept
{
    if (lambda <= zero)
        return 0;
    boost::random::poisson_distribution<uint> poissonDistr(lambda);
    return poissonDistr(g);
}


// Returns a binomially distributed deviate.
template <typename realT>
uint Boost<realT>::
binomial_number(
    const uint n,
    const realT p
)  noexcept
{
    boost::random::binomial_distribution<> bin_d(n, p);
    return bin_d(g);
}


// returns a multinomially distributed deviate
template <typename realT>
std::vector<uint> Boost<realT>::
multinomial_number(
    const uint n,
    const uint k
)  noexcept
{
    std::vector<realT> p (k - 1, one / k);
    return multinomial_number(n, p);
}


// Distributes n into p.size()+1 slots with probabilities
// p[0], p[1], ..., p.back(), 1 - sum(p)
template <typename realT>
std::vector<uint> Boost<realT>::
multinomial_number(
    const uint n,
    std::vector<realT> p   // by value
) noexcept
{
    std::vector<uint> res(p.size()+1);
    uint rem {n};
    for (szt i=0; i<p.size(); i++) {
        boost::binomial_distribution<> bin_d(rem, p[i]);
        res[i] = bin_d(g);
        rem -= res[i];
        for (szt j=i+1; j<p.size(); j++)
            p[j] /= p[i];
    }
    res.back() = rem;

    return res;
}


// Weibull distribution
template <typename realT> constexpr
realT Boost<realT>::
weibull_number(
    const realT lambda,
    const realT k
)  noexcept
{
    return lambda * std::pow((-std::log(one - r01u())),
                             one / k);
}


// Logistic distribution.
template <typename realT> constexpr
realT Boost<realT>::
logistic_number(
    const realT mi,
    const realT s
)  noexcept
{
    const auto u = r01u();
    return mi + s * (std::log(u) - std::log(one - u));
}


// Returns a normally distributed deviate N(mi, sigma^2).
template <typename realT> constexpr
realT Boost<realT>::
gaussian_number(
    const realT mi,
    const realT sigma
) noexcept
{
    return sigma * normalDistr(g) + mi;
}


template <typename realT> constexpr
realT Boost<realT>::
gaussian_number_constrained(
    const realT mi,
    const realT sigma,
    const realT cmin,
    const realT cmax
) noexcept
{
    XASSERT(cmin < cmax, "Error in gaussianNumConstr: cmin >= cmax.");
    
    realT res;
    do res = gaussian_number(mi, sigma);
    while (res < cmin || res > cmax);
    
    return res;
}

}    // namespace utils::random

#endif // UTILS_RANDOM_WITH_BOOST_H
