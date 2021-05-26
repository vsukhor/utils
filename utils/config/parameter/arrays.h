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
 * \file arrays.h
 * \brief Parameters of type std::array.
 * \details Contains template partial specialization for classes encapsulating
 * confuguration file parameter of type std::array and the exception handler
 * class for out of range errors.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_ARRAYS_H
#define UTILS_CONFIG_PARAMETER_ARRAYS_H

#include <array>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "../../common/misc.h"
#include "../../common/msgr.h"
#include "../exceptions/arrays.h"
#include "base.h"

/// Namespace for config parameters.
namespace Utils::Config::Parameter {


/**
* \brief parameters of std arrays of continuous fundamental types.
* \details Partial template specialization for for parameters of std arrays of
* continuous fundamental types.
* \tparam T Parameter type: must be std::is_fundamental.
* \tparam W Length of the std::array.
*/
template <typename T, szt W>
class Par<std::array<T,W>, false,
          std::enable_if_t<std::is_fundamental<T>::value>>
    : public Base<T> {

    using Q = std::array<T,W>;

    using Base<T>::get_name;
    using Base<T>::isLoaded_;

    Q p_;        ///< The parameter value.
    
public:

    /**
    * \brief Constructor.
    * \param name Name of the parameter.
    */
    explicit Par(const std::string& name);

    /**
    * \brief Constructor.
    * \param name Name of the parameter.
    * \param fname Configuration file name.
    * \param range Acceptable range of parameter values.
    * \param msgr \a Msgr used for the output.
    * \see Msgr 
    */
    explicit Par(const std::string& name,
                  const std::string& fname,
                  const std::vector<Q>& range,
                  Msgr* msgr=nullptr);

    /**
    * \brief Check that the read in parameter value is within the range set by \p r.
    * \param r Acceptable range of parameter values.
    * \param msgr \a Msgr used for the output.
    * \see Msgr
    */
    void check_range(const std::vector<Q>& r,
                     Msgr* msgr=nullptr);

    /**
    * \brief Print the the parameter to std::cout and logfile.
    * \param msgr \a Msgr used for the output.
    * \see Msgr
    */
    void print(Msgr* msgr=nullptr) final;

    /**
    * \brief The parameter values.
    * \return Parameter values (the whole array).
    */
    Q operator()() const;

    /**
    * \brief Specific component of the parameter array.
    * \param i Index in the array.
    * \return Parameter value (the \p i -th component).
    */
    T operator[](const szt i) const;

private:
    
    /**
    * \brief Initialize the parameter from the config file.
    * \param value Value to search for.
    */
    void initialize(std::string value) final;
};    

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T, szt W>
Par<std::array<T,W>, false,
    std::enable_if_t<std::is_fundamental<T>::value>>::
Par( const std::string& name )
    : Base<T> {name}
{}

template <typename T, szt W>
Par<std::array<T,W>, false,
    std::enable_if_t<std::is_fundamental<T>::value>>::
Par( const std::string& name,
     const std::string& fname,
     const std::vector<Q>& range,
     Msgr* msgr )
    : Base<T> {name}
{
    this->load(fname);
    try {
        check_range(range, msgr);
    } catch (const Exceptions::ParOutOfRange<T,false>&) {
        exit(0);
    }
    print(msgr);
}

template <typename T, szt W>
void Par<std::array<T,W>, false,
         std::enable_if_t<std::is_fundamental<T>::value>>::
check_range( const std::vector<Q>& r,
             Msgr* msgr )
{
    for (szt i=0; i<W; i++)
        if (p_[i]<r[0][i] || p_[i]>r[1][i])
            throw Exceptions::ParOutOfRange<T,false> {
                get_name(), p_[i], std::array<T,2>{r[0][i], r[1][i]}, msgr
            };
}

template <typename T, szt W>
void Par<std::array<T,W>, false,
         std::enable_if_t<std::is_fundamental<T>::value>>::
print( Msgr* msgr )
{
    (msgr)
    ? msgr->print_array(get_name(), p_)
    : [&]() {
         std::cout << get_name() << " ";
         for (const auto o : p_)
             std::cout << o << " ";
        std::cout<< std::endl;
      }();
}

template <typename T, szt W>
std::array<T,W> Par<std::array<T,W>, false,
                    std::enable_if_t<std::is_fundamental<T>::value>>::
operator()() const
{
    XASSERT(true, get_name());
    return p_;
}

template <typename T, szt W>
T Par<std::array<T,W>, false,
      std::enable_if_t<std::is_fundamental<T>::value>>::
operator[]( const szt i ) const
{
    XASSERT(isLoaded_, get_name());
    XASSERT(i<W, get_name());
    return p_[i];
}

template <typename T, szt W>
void Par<std::array<T,W>, false,
         std::enable_if_t<std::is_fundamental<T>::value>>::
initialize( std::string value )
{
    const std::string emp {" "};
    const std::string tab {"\t"};
    szt i {};
    while (value.length() && i<W) {
        ulong e {value.find(emp)};
        if (e == std::string::npos) e = value.find(tab);
        if (e == std::string::npos) e = value.length();
        const std::string val {value.substr(0, e)};
        if (val.length() < 1)
            throw Common::Exceptions::Simple
                {"Error in config file: Number of elelments in " + get_name() +
                 " is " + Common::STR(p_.size()) + " which is insufficient",
                 nullptr};
        std::stringstream(val) >> p_[i];
        value.erase(0, e);
        while (!value.substr(0, 1).compare(emp) ||
               !value.substr(0, 1).compare(tab))
            value.erase(value.begin());
        i++;
    }
    if (i < W)
        throw Common::Exceptions::Simple {
            "Improper Config::" + get_name() +
            " initialization: Data size insufficient", nullptr};
}

}  // namespace Utils::Config::Parameter

#endif // UTILS_CONFIG_PARAMETER_ARRAYS_H
