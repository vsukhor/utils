#include <filesystem>

#include "gtest/gtest.h"

#include "../utils/config/parameter/arrays.h"
#include "../utils/config/parameter/scalars.h"
#include "../utils/config/parameter/vectors.h"
#include "../utils/config/exceptions/arrays.h"
#include "../utils/config/exceptions/scalars.h"
#include "../utils/config/exceptions/vectors.h"

namespace {

class ParameterTest
    : public testing::Test {

protected:

    ParameterTest()
        : file {std::filesystem::current_path() / "tests" /
                       "data" / "config_sample.txt" }
    {}

    const std::filesystem::directory_entry file;
};

TEST_F(ParameterTest, Long1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<long, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, Unsigned1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<unsigned, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, Int1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<int, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, Double1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<double, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, ArrayInt1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::array<int,3>, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, VectorInt1)
{
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::vector<int>, false>;

    const std::string good_name {"some_name"};
    const std::string bad_name {"some name"};

    Par p {good_name};
    EXPECT_STREQ(p.get_name().c_str(), good_name.c_str());
    ASSERT_EXIT(Par {bad_name}, testing::ExitedWithCode(EXIT_FAILURE),
        "Incorrect name of configuration parameter: " + bad_name);
}

TEST_F(ParameterTest, Int2)
{
    using T = int;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<T, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    Par p1 {"intPositive", file, utils::common::zerohuge<T>};
    ASSERT_EQ(p1(), 74);
    Par p2 {"badUInt", file, utils::common::mhugehuge<T>};
    ASSERT_EQ(p2(), -74);
    ASSERT_EXIT(Par("badUInt", file, utils::common::zerohuge<T>),
        testing::ExitedWithCode(EXIT_FAILURE), "");
    testing::internal::GetCapturedStdout();
}

TEST_F(ParameterTest, Double2)
{
    using T = double;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<T, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    Par p1 {"myValue", file, utils::common::zerohuge<T>};
    ASSERT_EQ(p1(), 30.2);
    Par p2 {"negativeReal", file, utils::common::mhugehuge<T>};
    ASSERT_EQ(p2(), -45.3e-6);
    ASSERT_EXIT(Par("negativeReal", file, utils::common::zerohuge<T>),
        testing::ExitedWithCode(EXIT_FAILURE), "");
    testing::internal::GetCapturedStdout();
}


TEST_F(ParameterTest, VectorDouble2)
{
    using T = double;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::vector<T>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {0.87, 33.4, 64.326, 50.22};
    Par p1 {"reaction_rates", 4};
    p1.load(file);
    for (int i=0; i<rr.size(); i++)
        ASSERT_DOUBLE_EQ(p1()[i], rr[i]);
    const std::vector<std::array<T,2>> range1(4, {{0., 2000.}});
    const std::vector<std::array<T,2>> range2(4, {{0., 2.}});
    EXPECT_NO_THROW(p1.check_range(range1, nullptr));
    EXPECT_ANY_THROW(p1.check_range(range2, nullptr));
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

TEST_F(ParameterTest, VectorInt2)
{
    using T = int;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::vector<T>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {2, 12};
    Par p1 {"dimensionality", 2};
    p1.load(file);
    for (int i=0; i<rr.size(); i++)
        ASSERT_DOUBLE_EQ(p1()[i], rr[i]);
    const std::vector<std::array<T,2>> range1(2, {{0, 2000}});
    const std::vector<std::array<T,2>> range2(2, {{0, 2}});
    EXPECT_NO_THROW(p1.check_range(range1, nullptr));
    EXPECT_ANY_THROW(p1.check_range(range2, nullptr));
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

TEST_F(ParameterTest, VectorDouble3)
{
    using T = double;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::vector<T>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {0.87, 33.4, 64.326, 50.22};
    const std::vector<std::array<T,2>> range1(4, {{0., 2000.}});
    Par p1 {"reaction_rates", file, range1, nullptr};
    for (int i=0; i<rr.size(); i++)
        ASSERT_DOUBLE_EQ(p1()[i], rr[i]);
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

TEST_F(ParameterTest, VectorInt3)
{
    using T = int;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::vector<T>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {2, 12};
    const std::vector<std::array<T,2>> range1(2, {{0, 2000}});
    Par p1 {"dimensionality", file, range1, nullptr};
    for (int i=0; i<rr.size(); i++)
        ASSERT_EQ(p1()[i], rr[i]);
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

TEST_F(ParameterTest, ArrayDouble2)
{
    using T = double;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::array<T,4>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {0.87, 33.4, 64.326, 50.22};
    const std::vector<std::array<T,4>> range1{
        {{0., 0., 0., 0.}}, {{2000., 2000., 2000., 2000.}}
    };
    EXPECT_NO_THROW(Par("reaction_rates", file, range1, nullptr));
    Par p1 {"reaction_rates", file, range1, nullptr};
    for (int i=0; i<rr.size(); i++)
        ASSERT_DOUBLE_EQ(p1()[i], rr[i]);
    const std::vector<std::array<T,4>> range2{
        {{0., 0., 0., 0.}}, {{2., 2., 2., 2.}}
    };
    EXPECT_ANY_THROW(p1.check_range(range2, nullptr));
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

TEST_F(ParameterTest, ArrayInt2)
{
    using T = int;
    namespace nsp = utils::config::parameter;
    using Par = nsp::Par<std::array<T,2>, false>;

    testing::internal::CaptureStdout(); // avoide polluting stdout
    testing::internal::CaptureStderr(); // avoide polluting stderr
    const auto rr = std::vector {2, 12};
    const std::vector<std::array<T,2>> range1{
        {{0, 0}}, {{2000, 2000}}
    };
    Par p1 {"dimensionality", file, range1, nullptr};
    for (int i=0; i<rr.size(); i++)
        ASSERT_EQ(p1()[i], rr[i]);
    const std::vector<std::array<T,2>> range2{{{0, 0}}, {{2, 2}}};
    EXPECT_NO_THROW(p1.check_range(range1, nullptr));
    EXPECT_ANY_THROW(p1.check_range(range2, nullptr));
    testing::internal::GetCapturedStdout();
    testing::internal::GetCapturedStderr();
}

}  // namespace
