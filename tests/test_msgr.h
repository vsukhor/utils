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
        : scratch_dir {std::filesystem::current_path()/"tests"/"scratch_msgr"}
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

}  // namespace
