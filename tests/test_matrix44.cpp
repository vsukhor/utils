/* =============================================================================

 Copyright (C) 2009-2023 Valerii Sukhorukov. All Rights Reserved.

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

//#include <array>
//#include <string>

#include "gtest/gtest.h"

#include "../utils/matrices/matrix44.h"

namespace {

using namespace utils;

using M = matrices::Matrix<double, 4>;

constexpr auto s {M::order};

/// Random 2D array 1:
constexpr double _m1[s][s]
    {{-4212.7104166,   9347.78011337, -9999.68122136, -5005.97324314},
     {  135.32679396, -6906.35565029, -1787.6950769,   2020.88264006},
     {-1492.86349775, -4806.09238766,  7926.19697613,  3016.20654532},
     { 4940.28243447, -7752.10676624,  8223.92409579,  8103.09200627}};

/// Random 2D array 2:
constexpr double _m2[s][s]
    {{-1074.44510427, -8246.7596346,   4974.83328852, -4532.36902063},
     { 2394.59919692, -4376.90425763, -1911.89314074, -8498.78679104}, 
     {-8887.97024634,  7495.90053062, -8975.63509124, -8515.07752035}, 
     { 9062.80881424, -9953.25275575, -3358.41522101, -6247.51288577}};

M m1 {_m1};
M m2 {_m2};


TEST(TestMatrix4, Constructor)
{
    M m {_m1};
    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++)
            ASSERT_DOUBLE_EQ(m(i,j), _m1[i][j]);
    
    M mm {m};
    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++)
            ASSERT_DOUBLE_EQ(mm(i,j), m(i,j));
}


TEST(TestMatrix4, OpEqual)
{
    M m {_m1};
    ASSERT_TRUE(m == m1);

    m(1,2) = m2(1,2);
    ASSERT_FALSE(m == m1);
}


TEST(TestMatrix4, OpSum)
{
    constexpr double _m[s][s]
        {{ -5287.15552086,   1101.02047877,  -5024.84793284,  -9538.34226377}, 
         {  2529.92599089, -11283.25990792,  -3699.58821764,  -6477.90415098}, 
         {-10380.83374408,   2689.80814296,  -1049.43811511,  -5498.87097503}, 
         { 14003.09124871, -17705.35952198,   4865.50887478,   1855.57912049}};
    M m {_m};

    M mm {m1 + m2};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++)
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix4, OpDif)
{
    constexpr double _m[s][s]
        {{ -3138.26531233,  17594.53974797, -14974.51450989,   -473.60422251}, 
         { -2259.27240296,  -2529.45139266,    124.19806385,  10519.6694311 }, 
         {  7395.10674859, -12301.99291828,  16901.83206737,  11531.28406567}, 
         { -4122.52637977,   2201.14598951,  11582.33931679,  14350.60489204}};
    M m {_m};

    M mm {m1 - m2};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix4, Scale)
{
    constexpr double _m[s][s]
        {{-421.27104166,  934.77801134, -999.96812214, -500.59732431}, 
         {  13.5326794,  -690.63556503, -178.76950769,  202.08826401}, 
         {-149.28634977, -480.60923877,  792.61969761,  301.62065453}, 
         { 494.02824345, -775.21067662,  822.39240958,  810.30920063}};
    M m {_m};

    M mm {m1.scale(0.1)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix4, OpEmeMul)
{
    constexpr double _m[s][s]
        {{  4526326.08280666, -77088895.71206717, -49746747.01466505,   22688918.04529586}, 
         {   324053.43214817,  30228457.4504807,    3417881.95525598,  -17175050.68759913}, 
         { 13268526.34980537, -36025990.47888255, -71142651.71901411,  -25683232.55074631}, 
         { 44772835.19194344,  77158678.03389831, -27619351.85970971,  -50624171.7237784 }};
    M m {_m};

    M mm {m1.ele_mul(m2)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-4);
}


TEST(TestMatrix4, OpMatMul)
{
    constexpr double _m[s][s]
        {{ 70419203.57189474, -31304027.1576101,   67736137.71202886,   56071711.15298997}, 
         { 17520498.75076655,  -4402290.27317985,  23136177.88112348,   60679165.03441618}, 
         {-53017164.72135751,  62740010.67743962, -79510337.64457306,  -38723808.18134432}, 
         {-23528469.926797,   -25817498.50138395, -61630207.97049666,  -77159203.3923422 }};
    M m {_m};

    M mm {m1.mat_mul(m2)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-3);
}


TEST(TestMatrix4, Transpose)
{
    constexpr double _m1[s][s]
        {{-4212.7104166,    135.32679396, -1492.86349775,  4940.28243447}, 
         { 9347.78011337, -6906.35565029, -4806.09238766, -7752.10676624}, 
         {-9999.68122136, -1787.6950769,   7926.19697613,  8223.92409579}, 
         {-5005.97324314,  2020.88264006,  3016.20654532,  8103.09200627}};
    M tm {_m1};

    M tm1 {m1.t()};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_DOUBLE_EQ(tm(i,j), tm1(i,j));
}


TEST(TestMatrix4, Det)
{
    M mm1 {m1.scale(0.001)};
    constexpr double d_mm1 {1542.89433648472};
    const double det_mm1 {mm1.det()};
    ASSERT_NEAR(d_mm1, det_mm1, 1.e-3);

    M mm2 {m2.scale(0.001)};
    constexpr double d_mm2 {-4565.557085516252};
    const double det_mm2 {mm2.det()};
    ASSERT_NEAR(d_mm2, det_mm2, 1.e-3);
}


TEST(TestMatrix4, Inv)
{
    constexpr double _k[s][s]
        {{-0.1657814,  -0.06685971, -0.2203784,  -0.00371153}, 
         { 0.09519296, -0.11572999,  0.00493435,  0.08583478}, 
         {-0.07595581, -0.09147459,  0.05657213, -0.04516875}, 
         { 0.26923151,  0.02288445,  0.08166491,  0.25363167}};
    M k {_k};

    M m1s {m1.scale(0.001)};
    M im1;
    m1s.invert(im1);

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(im1(i,j), k(i,j), 1.e-7);
}


TEST(TestMatrix4, Diag)
{
    std::array<double, s> _d
        {{-4212.7104166,  -6906.35565029,  7926.19697613,  8103.09200627}};
    
    std::array<double, s> d {m1.diag()};

    for (int i=0; i<s; i++)
        ASSERT_DOUBLE_EQ(_d[i], d[i]);
}


TEST(TestMatrix4, Trace)
{
    constexpr double tr_mm1 {4910.222915509999};
    const double tr_m1 {m1.trace()};
    ASSERT_NEAR(tr_mm1, tr_m1, 1.e-7);

    constexpr double tr_mm2 {-20674.497338909998};
    const double tr_m2 {m2.trace()};
    ASSERT_NEAR(tr_mm2, tr_m2, 1.e-7);
}


}  // namespace
