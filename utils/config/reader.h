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
 * \file reader.h
 * \brief Convenience class reading parameters from confuguration files.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_READER_H
#define UTILS_CONFIG_READER_H

#include <array>
#include <filesystem>
#include <vector>

#include "../common/misc.h"
#include "../common/msgr.h"
#include "parameter/arrays.h"
#include "parameter/scalars.h"
#include "parameter/vectors.h"

/// Configuration namespace.
namespace utils::config {

/**
* Convenience class reading parameters from confuguration files.
*/
class Reader {

    using Msgr = common::Msgr;
    using path = std::filesystem::path;

public:

    const path file;  ///< Name of the configuration file.

    /**
    * \brief Constructor creating the configuration file-specific instance.
    * \param file Name of the configuration file.
    * \param msgr Messanger used for outputing.
    */
    Reader( const path& file,
            Msgr* msgr
        )
        : file {check_name(file)}
        , msgr {msgr}
    {
        if (msgr != nullptr)
            msgr->print("\nReading config from: ", file);
    }
    
    /**
    * \brief Value of the parameter read in from the configuration file.
    * \tparam T Parameter type.
    * \param s Parameter name.
    * \param range Acceptable range of parameter values.
    * \return Value of the parameter read in.
    */
    template <typename T>
    auto operator()(const std::string& s,
                    const std::vector<T>& range ) const
    {
        return parameter::Par<T,true>(s, file, range, msgr)();
    }

    /**
    * \brief Value of the parameter read in from the configuration file.
    * \tparam T Parameter type.
    * \param s Parameter name.
    * \param range Acceptable range of parameter values.
    * \return Value of the parameter read in.
    */
    template <typename T>
    auto operator()(const std::string& s,
                    const std::array<T,2>& range ) const
    {
        return parameter::Par<T,false>(s, file, range, msgr)();
    }

    /**
    * \brief Value of the parameter read in from the configuration file.
    * \tparam T Parameter type.
    * \param s Parameter name.
    * \param range Acceptable range of parameter values.
    * \return Value of the parameter read in.
    */
    template <typename T, auto N>
    auto operator()(const std::string& s,
                    const std::vector<std::array<T,N>>& range) const
    {
        return parameter::Par<std::array<T,N>,false>(s, file, range, msgr)();
    }

    /**
    * \brief Checks if file with name \p fname exists.
    * \param f Expected file.
    * \return Name of the confuguration file if it is found.
    */
    static auto check_name(const path& f) -> path
    {
        if (std::filesystem::exists(f) &&
              std::filesystem::is_regular_file(f))
            return f;

        throw common::exceptions::Simple(
            "Error: file provided '" + f.string() +
            "' is not a valid config file.", nullptr);
    }

    /**
    * \brief Copies \p compartment - specific config. file to directory \p path.
    * \note The directory is expected to exist.
    * \param path Directory name to which file should be copied.
    * \param signature Case-specific signature present in the file name.
    * \param compartment Name of the compartment specified in the configuration.
    */
    void copy( const std::filesystem::path& path,
               const std::string& signature,
               const std::string& compartment ) const
    {
        const std::filesystem::path cfgCopy {
            path / (std::string("cfgCopy_")+compartment+signature+".txt")
        };
        msgr->print("Copying "+compartment+" config to ", cfgCopy);
        if (std::filesystem::exists(cfgCopy)) {
            try {
                std::filesystem::remove(cfgCopy);
            }
            catch (const std::filesystem::filesystem_error& e) {
                msgr->print(e.what());
            }
        }
        try {
            std::filesystem::copy_file(file, cfgCopy);
        }
        catch (const std::filesystem::filesystem_error& e) {
                msgr->print(e.what());
        }
    }

private:

    Msgr* msgr {};  ///< \a Msgr output message processor.
};

}  // namespace utils::config

#endif // UTILS_CONFIG_READER_H
