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
 * \file vectors.h
 * \brief Parameters of type std::vector.
 * \details Contains template partial specialization for classes encapsulating
 * confuguration file parameter of type std::vector.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_VECTORS_H
#define UTILS_CONFIG_PARAMETER_VECTORS_H

#include <filesystem>
#include <sstream>
#include <type_traits>
#include <vector>

#include "../../common/misc.h"
#include "../../common/msgr.h"
#include "base.h"

namespace utils::config::parameter {

// specialization for vectors of fundamental types xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
* \brief Parameters of std vector of continuous fundamental types.
* \details Partial template specialization for for parameters of std vector of
* continuous fundamental types.
* \tparam T Parameter type: must be std::is_fundamental.
* \tparam isDiscrete Specifies if the vector components accept discrete
* or continous values.
*/
template <typename T, bool isDiscrete>
class Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_fundamental<T>::value>>
    : public Base<T>
{
    using Q = std::vector<T>;

    using Base<T>::get_name;
    using Base<T>::isLoaded_;

    Q p_;                ///< The parameter value.
    szt expectedSize_ {common::huge<szt>};    ///< Expected size of the vector.
    
public:    

    /**
    * \brief Constructor.
    * \param name Name of the parameter.
    * \param expectedSize Expected size of the parameter vector.
    */
    explicit Par(const std::string& name,
                 szt expectedSize);
                 
    /**
    * \brief Constructor.
    * \param name Name of the parameter.
    * \param file Configuration file.
    * \param range Acceptable range of parameter values.
    * \param msgr \a Msgr used for the output.
    * \see Msgr
    */
    explicit Par(const std::string& name,
                  const std::filesystem::directory_entry& file,
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
    * \brief Read a parameter without instantiating.
    * \details Static function for reading a parameter without instantiating
    * this class object.
    * \param name Name of the parameter.
    * \param file Configuration file.
    * \param msgr \a Msgr used for the output.
    * \return Parameter values (the whole vector).
    */
    static auto readin(
        const std::string& name,
        const std::filesystem::directory_entry& file,
        Msgr* msgr=nullptr
    );
    
    /**
    * \brief Print the the parameter to std::cout and logfile.
    * \param msgr \a Msgr used for the output.
    * \see Msgr
    */
    void print(Msgr* msgr=nullptr) final;
    
    /**
    * \brief The parameter values.
    * \return Parameter values (the whole vector).
    */
    Q operator()() const;
    
    /**
    * \brief Specific component of the parameter vector.
    * \param i Index in the vaector.
    * \return  Parameter value (the \p i -th component).
    */
    T operator[](szt i) const;

private:
    
    /**
    * \brief Initialize the parameter from the config file.
    * \param value Value to search for.
    */
    void initialize(std::string value) final;
};    


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T, bool isDiscrete>
Par<std::vector<T>,
    isDiscrete,
    std::enable_if_t<std::is_fundamental<T>::value>>::
Par( const std::string& name,
     const szt expectedSize )
    : Base<T> {name}
    , expectedSize_ {expectedSize}
{}


template <typename T, bool isDiscrete>
Par<std::vector<T>,
    isDiscrete,
    std::enable_if_t<std::is_fundamental<T>::value>>::
Par( const std::string& name,
     const std::filesystem::directory_entry& file,
     const std::vector<Q>& range,
     Msgr* msgr )
    : Base<Q> {name}
{
    this->load(file);
    check_range(range, msgr);
    print(msgr);
}


template <typename T, bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_fundamental<T>::value>>::
check_range(
    const std::vector<Q>& r,
    Msgr* msgr
)
{
    if (!r.size()) return;   // use this case to omit string checkups

    XASSERT(!isDiscrete || r.size()==2,
            "size of r must be 2 for continuous parameters");
    if constexpr (isDiscrete) {
        if (std::find(r.begin(), r.end(), p_) == r.end())
            throw exceptions::ParOutOfRange<Q,isDiscrete>
                {get_name(), p_, r, msgr};
    }
    else {
        if (p_<r[0] || p_>r[1])
            throw exceptions::ParOutOfRange<Q,isDiscrete>
                {get_name(), p_, r, msgr};
    }
}


template <typename T, bool isDiscrete>
auto Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_fundamental<T>::value>>::
readin(
    const std::string& name,
    const std::filesystem::directory_entry& file,
    Msgr* msgr
)
{
    return Par<Q,isDiscrete> {name, file, msgr}();
}


template <typename T, bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_fundamental<T>::value>>::
print( Msgr* msgr )
{
    (msgr != nullptr)
    ? msgr->print(get_name(), p_, 1)
    : std::cout << get_name() << " " << p_ << std::endl;
}


template <typename T, bool isDiscrete>
std::vector<T> Par<std::vector<T>,
                   isDiscrete,
                   std::enable_if_t<std::is_fundamental<T>::value>>::
operator()() const
{
    XASSERT(true, get_name());
    return p_;
}


template <typename T, bool isDiscrete>
T Par<std::vector<T>,
      isDiscrete,
      std::enable_if_t<std::is_fundamental<T>::value>>::
operator[]( const szt i ) const
{
    XASSERT(isLoaded_, get_name());
    XASSERT(i<p_.size(), get_name());
    return p_[i];
}


template <typename T, bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_fundamental<T>::value>>::
initialize( std::string value )
{
    const std::string emp {" "};
    const std::string tab {"\t"};
    const std::string exceptMessage
            {"Improper Config::" + get_name()+
             " initialization: Excessive data size"};

    while (value.length()) {
        ulong e {value.find(emp)};
        if (e == std::string::npos) e = value.find(tab);
        if (e == std::string::npos) e = value.length();
        const std::string val {value.substr(0, e)};
        if (val.length() < 1)
            throw common::exceptions::Simple
                    {"Error in config file: Number of elelments in " + get_name() +
                     " is " + STR(p_.size()) + " which is insufficient"};
        T tmp;
        std::stringstream(val) >> tmp;
        p_.push_back(tmp);
        value.erase(0, e);
        while (value.substr(0, 1) != emp ||
               value.substr(0, 1) != tab)
            value.erase(value.begin());
        if (p_.size() > expectedSize_)
            throw common::exceptions::Simple {exceptMessage, nullptr};
    }
    if (p_.size() != expectedSize_)
        throw common::exceptions::Simple {exceptMessage, nullptr};
}

}    // namespace utils::config::parameter

#endif // UTILS_CONFIG_PARAMETER_VECTORS_H
