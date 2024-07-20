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

#include "../utils/matrices/matrix33.h"

namespace {

using namespace utils;

/// Random 2D array 1:
constexpr double _m1[3][3]
    {{-4212.7104166,   9347.78011337, -9999.68122136},
     {  135.32679396, -6906.35565029, -1787.6950769},
     {-1492.86349775, -4806.09238766,  7926.19697613}};

/// Random 2D array 2:
constexpr double _m2[3][3]
    {{-1074.44510427, -8246.7596346,   4974.83328852},
     { 2394.59919692, -4376.90425763, -1911.89314074}, 
     {-8887.97024634,  7495.90053062, -8975.63509124}};

using M = matrices::Matrix<double, 3>;

constexpr auto s {M::order};

M m1 {_m1};
M m2 {_m2};


TEST(TestMatrix3, Constructor)
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


TEST(TestMatrix3, OpEqual)
{
    M m {_m1};
    ASSERT_TRUE(m == m1);

    m(1,2) = m2(1,2);
    ASSERT_FALSE(m == m1);
}


TEST(TestMatrix3, OpSum)
{
    constexpr double _m[s][s]
        {{ -5287.15552087,   1101.02047877,  -5024.84793284}, 
         {  2529.92599088, -11283.25990792,  -3699.58821764}, 
         {-10380.83374409,   2689.80814296,  -1049.43811511}};
    M m {_m};

    M mm {m1 + m2};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++)
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix3, OpDif)
{
    constexpr double _m[s][s]
        {{ -3138.26531233,  17594.53974797, -14974.51450989}, 
         { -2259.27240296,  -2529.45139266,    124.19806385 }, 
         {  7395.10674859, -12301.99291828,  16901.83206737}};
    M m {_m};

    M mm {m1 - m2};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix3, Scale)
{
    constexpr double _m[s][s]
        {{-421.27104166,  934.77801134, -999.96812214}, 
         {  13.5326794,  -690.63556503, -178.76950769}, 
         {-149.28634977, -480.60923877,  792.61969761}};
    M m {_m};

    M mm {m1.scale(0.1)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-7);
}


TEST(TestMatrix3, OpEmeMul)
{
    constexpr double _m[s][s]
        {{  4526326.08280666, -77088895.71206717, -49746747.01466505}, 
         {   324053.43214817,  30228457.4504807,    3417881.95525598}, 
         { 13268526.34980537, -36025990.47888255, -71142651.71901411}};
    M m {_m};

    M mm {m1.ele_mul(m2)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-4);
}


TEST(TestMatrix3, OpMatMul)
{
    constexpr double _m[s][s]
        {{ 1.15787382e+08, -8.11297441e+07,  5.09240010e+07}, 
         {-7.94374252e+05,  1.57120654e+07,  2.99231409e+07}, 
         {-8.03524680e+07,  9.27610768e+07, -6.93806637e+07}};
    M m {_m};

    M mm {m1.mat_mul(m2)};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(m(i,j), mm(i,j), 1.e-1);
}


TEST(TestMatrix3, Transpose)
{
    constexpr double _m1[s][s]
        {{-4212.7104166,    135.32679396, -1492.86349775}, 
         { 9347.78011337, -6906.35565029, -4806.09238766}, 
         {-9999.68122136, -1787.6950769,   7926.19697613}};
    M tm {_m1};

    M tm1 {m1.t()};

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_DOUBLE_EQ(tm(i,j), tm1(i,j));
}


TEST(TestMatrix3, Det)
{
    M mm1 {m1.scale(0.001)};
    constexpr double d_mm1 {391.3268651579034};
    const double det_mm1 {mm1.det()};
    ASSERT_NEAR(d_mm1, det_mm1, 1.e-3);

    M mm2 {m2.scale(0.001)};
    constexpr double d_mm2 {-479.2256094785707};
    const double det_mm2 {mm2.det()};
    ASSERT_NEAR(d_mm2, det_mm2, 1.e-3);
}


TEST(TestMatrix3, Inv)
{
    constexpr double _k[s][s]
        {{-0.16184159, -0.06652483, -0.21918336}, 
         { 0.00407884, -0.12347461, -0.02270293}, 
         {-0.02800891, -0.08739915,  0.07111567}};
    M k {_k};

    M m1s {m1.scale(0.001)};
    M im1;
    m1s.invert(im1);

    for (int i=0; i<s; i++)
        for (int j=0; j<s; j++) 
            ASSERT_NEAR(im1(i,j), k(i,j), 1.e-7);
}


TEST(TestMatrix3, Diag)
{
    std::array<double, s> _d
        {{-4212.7104166,  -6906.35565029,  7926.19697613}};
    
    std::array<double, s> d {m1.diag()};

    for (int i=0; i<s; i++)
        ASSERT_DOUBLE_EQ(_d[i], d[i]);
}


TEST(TestMatrix3, Trace)
{
    constexpr double tr_mm1 {-3192.8690907600003};
    const double tr_m1 {m1.trace()};
    ASSERT_NEAR(tr_mm1, tr_m1, 1.e-7);

    constexpr double tr_mm2 {-14426.98445314};
    const double tr_m2 {m2.trace()};
    ASSERT_NEAR(tr_mm2, tr_m2, 1.e-7);
}


}  // namespace
