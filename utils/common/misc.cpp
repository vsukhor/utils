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

============================================================================== */

/**
* \file misc.h
* \brief A loose collection of functions of common use.
* \author Valerii Sukhorukov
*/

#include <filesystem>

#include "misc.h"

/// General stuff.
namespace utils::common {

/*
bool file_exists( const std::string& name )
{
    if (auto file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } 
    else return false;
}
bool fileExists( const std::string& name )
{
  class stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

bool directory_exists( const std::string& pathstrconst )
// checks if a given directory exists.
{
    std::string pathstr {pathstrconst};
    char lascch {pathstr.back()};
    if (lascch == SLASH.c_str()[0])
        pathstr.pop_back();
    const char* path {pathstr.c_str()};
    class stat info;
    if (stat(path, &info) != 0)      return false;
    else if (info.st_mode & S_IFDIR) return true;
    else                             return false;
}

void check_directory( const std::string& s )
{
    namespace fs = std::filesystem;
    const fs::path path {s};
    if (!fs::exists(path))
        XASSERT(false, "Error: cannot create directory"+path.string());
    return;
//        if(!fs::create_directory(path))
//            XASSERT(false, "Error: cannot create directory"+path.string());
}
*/
/*
void copy_text_file(
    const std::string& from,
    const std::string& to
)
{
    std::ifstream ifs {from};
    XASSERT(ifs.is_open(), "Unable to open file for reading: "+from);
    std::ofstream ofs {to};
    XASSERT(ofs.is_open(), "Unable to open file for writing: "+to);

    for (std::string line; std::getline(ifs, line);)
        ofs << line << std::endl;
}
*/
//template float avg<float>( std::vector<float> const& );
//template double avg<double>( std::vector<double> const& );

}  // namespace utils::common
