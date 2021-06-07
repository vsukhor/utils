#include <string>

#include "gtest/gtest.h"

#include "../utils/common/misc.h"

namespace misc_test {

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

}
