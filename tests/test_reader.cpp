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

#include <filesystem>

#include "gtest/gtest.h"

#include "../utils/config/reader.h"
#include "../utils/constants.h"

namespace {

class ReaderTest
    : public testing::Test {

protected:

    using Reader = utils::config::Reader;
    static constexpr const int PRINT_PRECISION {6};

    const std::filesystem::path file;

    ReaderTest()
        : file {std::filesystem::current_path() / "tests" /
                       "data" / "config_sample.txt" }
    {}
};

TEST_F(ReaderTest, CheckName)
{
    EXPECT_NO_THROW(Reader::check_name(file));
    const std::filesystem::path fakefile {
        std::filesystem::current_path() / "tests" / "data" / "fake.txt" };
    testing::internal::CaptureStderr();
    EXPECT_ANY_THROW(Reader::check_name(fakefile));
    testing::internal::GetCapturedStderr();
}

TEST_F(ReaderTest, Constructor)
{
    EXPECT_NO_THROW(Reader(file, nullptr));
    Reader r {file, nullptr};
    ASSERT_STREQ(file.string().c_str(), r.file.string().c_str());
    const std::filesystem::path fakefile {
        std::filesystem::current_path() / "tests" / "data" / "fake.txt" };
    testing::internal::CaptureStderr();
    EXPECT_ANY_THROW(Reader(fakefile, nullptr));
    testing::internal::GetCapturedStderr();
}

TEST_F(ReaderTest, FunctScalar)
{
    Reader read {file, nullptr};
    testing::internal::CaptureStdout();
    ASSERT_DOUBLE_EQ(read("myValue", std::array<double,2>{{0., 100.}}), 30.2);
    ASSERT_DOUBLE_EQ(read("negativeReal", std::array<double,2>{{-10., 10.}}), -45.3e-6);
    ASSERT_DOUBLE_EQ(read("edge_length", std::array<double,2>{{0., 100'000.}}), 3.e+4);
    ASSERT_DOUBLE_EQ(read("badUInt", std::array<int,2>{{-100, 0}}), -74);
    testing::internal::GetCapturedStdout();
    ASSERT_EXIT(Reader(file, nullptr)("badUInt", std::array<int,2>{{0, 100}}),
                testing::ExitedWithCode(EXIT_FAILURE), "");
    ASSERT_EXIT(Reader(file, nullptr)("my Value", std::array<double,2>{{0., 100.}}),
                testing::ExitedWithCode(EXIT_FAILURE), "");
}

TEST_F(ReaderTest, FunctArrayDouble)
{
    const auto rr = std::vector {0.87, 33.4, 64.326, 50.22};
    const std::vector<std::array<double,4>> range1{
        {{0., 0., 0., 0.}}, {{2000., 2000., 2000., 2000.}}
    };
    Reader read {file, nullptr};
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(read("reaction_rates", range1));
    const auto p1 {read("reaction_rates", range1)};
    testing::internal::GetCapturedStdout();
    for (int i=0; i<rr.size(); i++)
        ASSERT_DOUBLE_EQ(p1[i], rr[i]);
}

TEST_F(ReaderTest, FunctArrayInt)
{
    const auto rr = std::vector {2, 12};
    const std::vector<std::array<int,2>> range1{
        {{0, 0}}, {{2000, 2000}}
    };
    Reader read {file, nullptr};
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(read("dimensionality", range1));
    const auto p1 {read("dimensionality", range1)};
    testing::internal::GetCapturedStdout();
    for (int i=0; i<rr.size(); i++)
        ASSERT_EQ(p1[i], rr[i]);
}

}  // namespace
