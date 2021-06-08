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

#include <filesystem>
#include <fstream>
#include <ostream>

#include "gtest/gtest.h"

#include "../utils/common/msgr.h"

namespace {

class MsgrTest
    : public testing::Test {

protected:

    using Msgr = utils::common::Msgr;
    static constexpr const int PRINT_PRECISION {6};

    MsgrTest()
        : scratch_dir {std::filesystem::current_path() / "tests" / "scratch"}
        , scratch_dir_created {}
        , file {scratch_dir / "msgr_test.txt"}
    {}

    void SetUp() override
    {
        if (!std::filesystem::exists(scratch_dir)) {
            scratch_dir_created = true;
            std::filesystem::create_directory(scratch_dir);
        }
        fstr = new Msgr::logstream {file};
    }

    void TearDown() override
    {
        if (fstr) {
            if (fstr->is_open()) fstr->close();
            delete fstr;
        }
        if (std::filesystem::exists(file))
            std::filesystem::remove(file);
        if (scratch_dir_created)
            std::filesystem::remove(scratch_dir);
    }

    const std::filesystem::path scratch_dir;
    bool scratch_dir_created;
    const std::filesystem::path file;
    Msgr::logstream* fstr {};
};

TEST_F(MsgrTest, ConstructorDefault)
{
    Msgr msgr {};

    ASSERT_EQ(msgr.so, nullptr);
    ASSERT_EQ(msgr.so, nullptr);
}

TEST_F(MsgrTest, Constructor)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};

    ASSERT_EQ(msgr.so, &std::cout);
    ASSERT_EQ(msgr.sl, fstr);
}

TEST_F(MsgrTest, SetFormats)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    int new_precision = 10;
    msgr.set_formats(new_precision);
    ASSERT_EQ(msgr.so->precision(), new_precision);
    ASSERT_EQ(msgr.sl->precision(), new_precision);

}

TEST_F(MsgrTest, PrintFalse)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"Hello world!"};
    std::string ss {"Hello world! "};

    testing::internal::CaptureStdout();
    msgr.print<false>(s0);
    std::string t1 {testing::internal::GetCapturedStdout()};
    ASSERT_STREQ(ss.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(ss.c_str(), t2.c_str());
}

TEST_F(MsgrTest, PrintTrue)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"Hello world!"};
    std::string ss {"Hello world! \n"};

    testing::internal::CaptureStdout();
    msgr.print<true>(s0);
    std::string t1 = testing::internal::GetCapturedStdout();
    ASSERT_STREQ(ss.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(ss.c_str(), (t2+"\n").c_str());
}

TEST_F(MsgrTest, PrintArrayULong)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"the name"};
    const std::array<unsigned long, 3> a0 {{34, 56, 78}};
    std::string exp = {s0 + "[]:   34 56 78  \n"};

    testing::internal::CaptureStdout();
    msgr.print_array<unsigned long,3>(s0, a0);
    std::string t1 = testing::internal::GetCapturedStdout();
    ASSERT_STREQ(exp.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(exp.c_str(), (t2+"\n").c_str());
}

TEST_F(MsgrTest, PrintArrayDouble)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"the name"};
    const std::array<double, 3> a0 {{34.74, 56.99, 78.62}};
    std::string exp = {s0 + "[]:   34.740000 56.990000 78.620000  \n"};

    testing::internal::CaptureStdout();
    msgr.print_array<double,3>(s0, a0);
    std::string t1 = testing::internal::GetCapturedStdout();
    ASSERT_STREQ(exp.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(exp.c_str(), (t2+"\n").c_str());
}

TEST_F(MsgrTest, PrintVectorULong)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"the name"};
    const std::vector<unsigned long> a0 {34, 56, 78};
    std::string exp = {s0 + "[]:   34 56 78 "};

    testing::internal::CaptureStdout();
    msgr.print_vector<unsigned long>(s0, a0);
    std::string t1 = testing::internal::GetCapturedStdout();
    ASSERT_STREQ(exp.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(exp.c_str(), t2.c_str());
}

TEST_F(MsgrTest, PrintVectorDouble)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s0 {"the name"};
    const std::vector<double> a0 {34.74, 56.99, 78.62};
    std::string exp = {s0 + "[]:   34.740000 56.990000 78.620000 "};

    testing::internal::CaptureStdout();
    msgr.print_vector<double>(s0, a0);
    std::string t1 = testing::internal::GetCapturedStdout();
    ASSERT_STREQ(exp.c_str(), t1.c_str());

    std::ifstream ifs {file};
    ASSERT_FALSE(ifs.fail());
    std::string t2 {};
    std::string line {};
    while (std::getline(ifs, line))
        t2 += line;
    ASSERT_STREQ(exp.c_str(), t2.c_str());
}

TEST_F(MsgrTest, Exit1)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s {"QWERTY"};

    testing::internal::CaptureStdout();
    ASSERT_EXIT(msgr.exit(s), testing::ExitedWithCode(EXIT_FAILURE), "");
    testing::internal::GetCapturedStdout();
}

TEST_F(MsgrTest, Exit2)
{
    Msgr msgr {&std::cout, fstr, PRINT_PRECISION};
    std::string s {"QWERTY"};

    testing::internal::CaptureStdout();
    ASSERT_EXIT(msgr.exit(s), testing::ExitedWithCode(EXIT_FAILURE), "");
    std::string t = testing::internal::GetCapturedStdout();
    ASSERT_STREQ((s+" \n").c_str(), t.c_str());
}

}  // namespace
