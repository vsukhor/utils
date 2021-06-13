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
 * \file base.h
 * \brief Contains the base class for configuration parameters.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_BASE_H
#define UTILS_CONFIG_PARAMETER_BASE_H


#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "../../common/misc.h"
#include "../../msgr.h"
#include "../../common/exceptions.h"

/// Configuration module
namespace utils::config {

using szt = common::szt;
using ulong = common::ulong;

namespace exceptions {

/**
 * Generic template for 'Parameter out of range' exception.
 * \tparam Q Parameter type.
 * \tparam isDiscrete Specifies if the parameter takes discrete values only.
 * \tparam Enabler SFINAE enabler.
 */
template <typename Q, bool isDiscrete, typename Enabler = void>
class ParOutOfRange
    : public std::exception
{};

}  // namespace exceptions


////////////////////////////////////////////////////////////////////////////////
namespace parameter {

/**
 * Base class for configuration parameters.
 * \tparam Q Parameter type.
 */
template <typename Q>
class Base {

public:

    using str = std::string;

    /**
     * Reads in parameters from input file stream \p ifs.
     * \param ifs Input file stream to load the parameters.
     */
    void load(std::ifstream& ifs);

    /**
     * Reads in parameters from file.
     * \param file Name of the configuration file to load the parameters.
     */
    void load(const std::filesystem::path& file);

    /**
     * Name of the parameter.
     * \return Name of the parameter.
     */
    auto get_name() const noexcept -> str;

protected:
    
    bool isLoaded_ {};  ///< Flag if the parameter is loaded.

    /**
     * Constructor.
     * \param name Name of the parameter.
     */
    explicit Base(str name);  // by pass-by-value + move

    // The rule of five is triggered by the virtual destructor,
    // the defaults suffice.
    Base(const Base&) = default;                ///< copy constructor
    Base& operator=(const Base&) = default;     ///< copy assignment
    Base(Base&&) = default;                     ///< move constructor
    Base& operator=(Base&&) noexcept = default; ///< move assignment
    virtual ~Base() noexcept = default;         ///< virtual destructor

    /**
     * Print the the parameter to std::cout and logfile.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    virtual void print(Msgr* msgr=nullptr) = 0;

    /**
     * Initialize the parameter from the config file.
     * \param value Value to search for.
     */
    virtual void initialize(str value) = 0;

    auto check_name(const str& s) const -> str;

private:

    const str name;   ///< parameter name

    /**
     * Finds the the parameter by \a name in the configuration file stream \p ifs.
     * \param ifs Input file stream to load the parameter.
     * \param[out] value str containig value(s) of the parameter searched.
     * \return Bool corresponding to the success/failure of the search.
     */
    bool detect_by_name(std::ifstream& ifs,
                        str& value) const;  // by reference
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Q>
Base<Q>::
Base( str name )
    : name {std::move(name)}
{}

// If the line contains a valid parname-value combination, returns true and the
// value, otherwise retruns false.
template <typename Q>
bool Base<Q>::
detect_by_name(
    std::ifstream& ifs,
    str& value
) const
{    
    const str emp {" "};
    const str tab {"\t"};
    
    str line;
    getline(ifs, line);
    
    ulong commentpos = line.find_first_of('#');
    if (commentpos != str::npos)
        line.erase(commentpos);
    
    if (!line.length())
        return false;
    while (!line.substr(line.length()-1, 1).compare(emp) ||
           !line.substr(line.length()-1, 1).compare(tab))
        line.erase(line.length()-1);
    if (!line.length())
        return false;
    
    int parnameend = -1;
    if (     line.find_first_of(emp) == str::npos &&
             line.find_first_of(tab) != str::npos)
        parnameend = static_cast<int>(line.find_first_of(tab));
    else if (line.find_first_of(emp) != str::npos &&
             line.find_first_of(tab) == str::npos)
        parnameend = static_cast<int>(line.find_first_of(emp));
    else if (line.find_first_of(emp) != str::npos &&
             line.find_first_of(tab) != str::npos)
        parnameend = std::min(static_cast<int>(line.find_first_of(emp)),
                              static_cast<int>(line.find_first_of(tab)));
    const auto parname = line.substr(0, static_cast<size_t>(parnameend));
    
    if (parname != name)
        return false; 
    
    value = line.substr(line.find_last_of("=") + 1);
    while (!value.substr(0, 1).compare(emp) ||
           !value.substr(0, 1).compare(tab))
        value.erase(value.begin());
    return true;
}

template <typename Q>
void Base<Q>::
load( const std::filesystem::path& file )
{
    str parname;
    str value;
    std::ifstream ifs {file};
    if (!ifs.is_open()) {
        throw common::exceptions::Simple
            {"Unable to open config file: " + file.string(), nullptr};
    }
    try {
        load(ifs);
    } catch (const common::exceptions::Simple&) {
        return;
    }
}

template <typename Q>
void Base<Q>::
load( std::ifstream& ifs )
{
    XASSERT(!isLoaded_, "Repeated load of " + name);
    
    ifs.clear();
    ifs.seekg(0, std::ios::beg);
    while (ifs.good()) {
        str value;
        if (!detect_by_name(ifs, value))
            continue;
        initialize(value);
        isLoaded_ = true;
        return;
    }
    throw common::exceptions::Simple
        {"Error: parameter not loaded: " + name, nullptr};
}

template <typename Q> inline
auto Base<Q>::
check_name( const str& s ) const -> str
{
    const auto ss {s.substr(0, s.find(" "))};
    if (s.size() != ss.size()) {
        std::cerr << "Incorrect name of configuration parameter: " + s;
        std::exit(EXIT_FAILURE);
    }
    return s;
}


template <typename Q> inline
auto Base<Q>::
get_name() const noexcept -> str
{
    return name;
}

// template for Par xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * \brief Generic template for parameter classes.
 * \tparam Q Parameter type.
 * \tparam isDiscrete Specifies if the parameter takes discrete values only.
*/
template <typename Q, bool isDiscrete, typename Enabler=void>
class Par : public Base<Q>
{};


}  // namespace parameter
}  // namespace utils::config

#endif // UTILS_CONFIG_PARAMETER_BASE_H
