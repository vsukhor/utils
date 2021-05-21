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
 * \file base.h
 * \brief Contains the base class for configuration parameters.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_BASE_H
#define UTILS_CONFIG_PARAMETER_BASE_H

#include <sstream>

#include "../../common/misc.h"
#include "../../common/msgr.h"
#include "../../common/exceptions.h"

/// Library-wide.
namespace Utils {
/// Configuration module
namespace Config {

namespace Exceptions {

using namespace Common;

/**
* \brief Generic template for 'Parameter out of range' exception.
* \tparam Q Parameter type.
* \tparam isDiscrete Specifies if the parameter takes discrete values only.
* \tparam Enabler SFINAE enabler.
*/
template <typename Q, bool isDiscrete, typename Enabler = void>
class ParOutOfRange
    : public std::exception
{};

}     // namespace Exceptions


////////////////////////////////////////////////////////////////////////////////
namespace Parameter {

using namespace Utils::Common;

/**
* \brief Base class for configuration parameters.
* \tparam Q Parameter type.
*/
template <typename Q>
class Base {

public:

    /**
    * \brief Reads in parameters from input file stream \p ifs.
    * \param ifs Input file stream to load the parameters.
    */
    void load(std::ifstream& ifs);

    /**
    * \brief Reads in parameters from file.
    * \param fname Name of the configuration file to load the parameters.
    */
    void load(const std::string& fname);

protected:
    
    bool isLoaded_ {};  ///< Flag if the parameter is loaded.

    /**
    * \brief Constructor.
    * \param name Name of the parameter.
    */
    explicit Base(const std::string& name);

    // The rule of five is triggered by the virtual destructor, the defaults suffice.
    Base(const Base&) = default;                ///< copy constructor
    Base& operator=(const Base&) = default;     ///< copy assignment
    Base(Base&&) = default;                     ///< move constructor
    Base& operator=(Base&&) = default;          ///< move assignment
    virtual ~Base() noexcept = default;         ///< virtual destructor

    /**
    * \brief Print the the parameter to std::cout and logfile.
    * \param msgr \a Msgr used for the output.
    * \see Msgr
    */
    virtual void print(Msgr* msgr=nullptr) = 0;

    /**
    * \brief Initialize the parameter from the config file.
    * \param value Value to search for.
    */
    virtual void initialize(std::string value) = 0;

    /**
    * \brief Name of the parameter.
    * \return Name of the parameter.
    */
    std::string get_name() const noexcept;

private:

    const std::string name;        ///< parameter name

    /**
    * \brief Finds the the parameter by \a name in the configuration file stream \p ifs.
    * \param ifs Input file stream to load the parameter.
    * \param[out] value std::string containig value(s) of the parameter searched.
    * \return Bool corresponding to the success/failure of the search.
    */
    bool detect_by_name(std::ifstream& ifs,
                        std::string& value) const;    // by reference
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename Q>
Base<Q>::
Base( const std::string& name )
    : name {name}
{}

// If the line contains a valid parname-value combination, returns true and the
// value, otherwise retruns false.
template <typename Q>
bool Base<Q>::
detect_by_name( std::ifstream& config,
                std::string& value ) const
{    
    const std::string emp {" "};
    const std::string tab {"\t"};
    
    std::string line;
    getline(config, line);
    
    ulong commentpos = line.find_first_of('#');
    if (commentpos != std::string::npos)
        line.erase(commentpos);
    
    if (!line.length())
        return false;
    while (!line.substr(line.length()-1, 1).compare(emp) ||
           !line.substr(line.length()-1, 1).compare(tab))
        line.erase(line.length()-1);
    if (!line.length())
        return false;
    
    int parnameend = -1;
    if (     line.find_first_of(emp) == std::string::npos &&
             line.find_first_of(tab) != std::string::npos)
        parnameend = static_cast<int>(line.find_first_of(tab));
    else if (line.find_first_of(emp) != std::string::npos &&
             line.find_first_of(tab) == std::string::npos)
        parnameend = static_cast<int>(line.find_first_of(emp));
    else if (line.find_first_of(emp) != std::string::npos &&
             line.find_first_of(tab) != std::string::npos)
        parnameend = std::min(static_cast<int>(line.find_first_of(emp)),
                              static_cast<int>(line.find_first_of(tab)));
    const auto parname = line.substr(0, static_cast<size_t>(parnameend));
    
    if (parname != name)
        return false; 
    
    value = line.substr(line.find_last_of("=")+1);
    while (!value.substr(0, 1).compare(emp) ||
           !value.substr(0, 1).compare(tab))
        value.erase(value.begin());
    return true;
}

template <typename Q>
void Base<Q>::
load( const std::string& fname )
{
    std::string parname, value;
    std::ifstream ifs {fname};
    if (!ifs.is_open()) {
        std::cout << "Unable to open config file: " << fname << std::endl;
        exit(0);
    }
    try {
        load(ifs);
    } catch (Utils::Common::Exceptions::Simple& e) {
        return;
    }
}

template <typename Q>
void Base<Q>::
load( std::ifstream& config )
{
    XASSERT(!isLoaded_, "Repeated load of "+name);
    
    config.clear();
    config.seekg(0, std::ios::beg);
    while (config.good()) {
        std::string value;
        if (!detect_by_name(config, value))
            continue;
        initialize(value);
        isLoaded_ = true;
        return;
    }
    throw Utils::Common::Exceptions::Simple
        {"Error: parameter not loaded: "+name};
}

template <typename Q>
std::string Base<Q>::
get_name() const noexcept
{
    return name;
}

// template for Par xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
* \brief Generic template for parameter classes.
* \tparam Q Parameter type.
* \tparam isDiscrete Specifies if the parameter takes discrete values only.
* \tparam Enabler SFINAE enabler
*/
template <typename Q, bool isDiscrete, typename Enabler=void>
class Par : public Base<Q>
{};


}    // namespace Parameter
}    // namespace Config
}   // namespace Utils

#endif // UTILS_CONFIG_PARAMETER_BASE_H
