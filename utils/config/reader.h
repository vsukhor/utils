/* ==============================================================================

 Copyright (C) 2009-2019, Valerii Sukhorukov, <vsukhorukov@yahoo.com>

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
 * \file reader.h
 * \brief Convenience class reading parameters from confuguration files.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_READER_H
#define UTILS_CONFIG_READER_H

#include "../common/misc.h"
#include "../common/msgr.h"
#include "parameter/vectors.h"
#include "parameter/arrays.h"
#include "parameter/scalars_strings.h"

/// Library-wide.
namespace Utils {
/// Configuration module
namespace Config{

using namespace Common;
using namespace Parameter;

/**
* Convenience class reading parameters from confuguration files.
*/
class Reader {

public:

	const std::string fname;		///< name of the configuration file.

	/**
	* \brief Constructor creating the configuration file-specific instance.
	* \param fname Name of the configuration file.
	* \param msgr Messanger used for outputing.
	*/
	Reader( const std::string& fname,
			Msgr* msgr
		)
		: fname {check_fname(fname)}
		, msgr {msgr}
	{
		if (msgr) msgr->print("\nReading config from: "+fname);
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
		return Par<T,true>(s, fname, range, msgr)();
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
		return Par<T,false>(s, fname, range, msgr)();
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
					const vecarr<T,N>& range) const
	{
		return Par<std::array<T,N>,false>(s, fname, range, msgr)();
	}

	/**
	* \brief Checks if file with name \p fname exists.
	* \param fname Expected name of the file.
	* \return Name of the confuguration file if it is found.
	*/
	static std::string check_fname( const std::string fname )
	{
		if (file_exists(fname))
			return fname;
		else
			exit(Utils::Common::Exceptions::simple("Error: no config file provided "+fname));
	}

	/**
	* \brief Copies \p compartment - specific configuration file to directory \p path.
	* \note The directory is expected to exist.
	* \param path Directory name to which file should be copied.
	* \param signature Case-specific signature present in the file name.
	* \param compartment Name of the compartment specified in the configuration.
	*/
	void copy( const std::string& path,
			   const std::string& signature,
			   const std::string& compartment ) const
	{
		const auto cfgCopy {path+"cfgCopy_"+compartment+signature+".txt"};
		msgr->print("Copying "+compartment+" config to "+cfgCopy);
		copy_text_file(fname, cfgCopy);
	}

	/**
	* \brief Extracts directory name.
	* \return Directory name.
	*/
	std::string path() const noexcept
	{
		auto found {fname.find_last_of("/\\")};
		XASSERT(found < fname.npos, "ConfigReader: fname cannot be splitted");
		return fname.substr(0, found)+SLASH;
	}

private:

	Msgr* msgr {};		///< \a Msgr output message processor.
};

}	// namespace Config
}   // namespace Utils

#endif // UTILS_CONFIG_READER_H
