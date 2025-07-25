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

#include <array>
#include <string>

#include "gtest/gtest.h"

#include "../utils/common/misc.h"

namespace {

TEST(TestMisc, Trim)
{
    const std::string s0 {" "};
    ASSERT_STREQ(utils::common::trim(s0).c_str(), std::string("").c_str());

    const std::string s1 {"  ab"};
    ASSERT_STREQ(utils::common::trim(s1).c_str(), std::string("ab").c_str());

    const std::string s2 {"ab   "};
    ASSERT_STREQ(utils::common::trim(s2).c_str(), std::string("ab").c_str());

    const std::string s3 {"    ab   "};
    ASSERT_STREQ(utils::common::trim(s3).c_str(), std::string("ab").c_str());
}


}  // namespace
