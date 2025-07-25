/* =============================================================================

 Copyright (C) 2009-2025 Valerii Sukhorukov. All Rights Reserved.

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

#include "../../constants.h"
#include "../../msgr.h"
#include "../exceptions/vectors.h"
#include "base.h"

namespace utils::config::parameter {

// specialization for vectors of arithmetic types xxxxxxxxxxxxxxxxxxxxxxxxxxxxx

/**
 * \brief Parameters of continuous arithmetic types packed into std vector.
 * \details Partial template specialization for for parameters of std vector of
 * continuous arithmetic types.
 * \tparam T Parameter type: must be std::is_arithmetic.
 * \tparam isDiscrete Specifies if the vector components accept discrete
 * or continous values.
 */
template<typename T,
         bool isDiscrete>
class Par<std::vector<T>,
          isDiscrete,
          std::enable_if_t<std::is_arithmetic_v<T>>>
    : public Base<T>
{
    using Q = std::vector<T>;

    using Base<T>::check_name;
    using Base<T>::isLoaded_;

    Q p_;  ///< The parameter value.

    szt expectedSize_ {undefined<szt>};   ///< Expected size of the vector.

public:

    using str = std::string;

    using Base<T>::get_name;
    /**
     * Constructor.
     * \param name Name of the parameter.
     */
    explicit Par(const str& name);

    /**
     * Constructor.
     * \param name Name of the parameter.
     * \param expectedSize Expected size of the parameter vector.
     */
    explicit Par(const str& name,
                 szt expectedSize);

    /**
     * Constructor.
     * \param name Name of the parameter.
     * \param file Configuration file.
     * \param range Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    explicit Par(const str& name,
                 const std::filesystem::path& file,
                 const std::vector<std::array<T, 2>>& range,
                 Msgr* msgr=nullptr);

    /**
     * \brief Checks that the parameter value is within the range set by \p r.
     * \param r Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    void check_range(const std::vector<std::array<T, 2>>& r,
                     Msgr* msgr=nullptr);

    /**
     * \brief Reads a parameter without instantiating.
     * \details Static function for reading a parameter without instantiating
     * this class object.
     * \param name Name of the parameter.
     * \param file Configuration file.
     * \param msgr \a Msgr used for the output.
     * \return Parameter values (the whole vector).
     */
    static auto readin(
        const str& name,
        const std::filesystem::path& file,
        Msgr* msgr=nullptr
    );

    /**
     * \brief Print the the parameter to std::cout and logfile.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    void print( Msgr* msgr=nullptr ) final;

    /**
     * \brief The parameter values.
     * \return Parameter values (the whole vector).
     */
    Q operator()() const;

    /**
     * \brief Accesses a specific component of the parameter vector.
     * \param i Index in the vector.
     * \return Parameter value (the \p i -th component).
     */
    T operator[]( szt i ) const;

private:

    /**
     * \brief Initialize the parameter from the config file.
     * \param value Value to search for.
     */
    void initialize( str value ) final;
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename T,
         bool isDiscrete>
Par<std::vector<T>,
    isDiscrete,
    std::enable_if_t<std::is_arithmetic_v<T>>>::
Par(
    const str& name
)
    : Base<T> {check_name(name)}
{}


template<typename T,
         bool isDiscrete>
Par<std::vector<T>,
    isDiscrete,
    std::enable_if_t<std::is_arithmetic_v<T>>>::
Par(
    const str& name,
    const szt expectedSize
)
    : Base<T> {check_name(name)}
    , expectedSize_ {expectedSize}
{}


template<typename T,
         bool isDiscrete>
Par<std::vector<T>,
    isDiscrete,
    std::enable_if_t<std::is_arithmetic_v<T>>>::
Par(
    const str& name,
    const std::filesystem::path& file,
    const std::vector<std::array<T, 2>>& range,
    Msgr* msgr
)
    : Base<T> {check_name(name)}
{
    this->load(file);
    try {
        check_range(range, msgr);
    }
    catch (const exceptions::ParOutOfRange<T, isDiscrete>&) {
        std::exit(EXIT_FAILURE);
    }
    print(msgr);
}


template<typename T,
         bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_arithmetic_v<T>>>::
check_range(
    const std::vector<std::array<T, 2>>& r,
    Msgr* msgr
)
{
    if (!r.size()) return;   // use this case to omit string checkups

    ASSERT(r.size() == p_.size(),
           "Size of range vector must be to the size of imported vector.");

    for (szt i=0; i<r.size(); i++)
        if (p_[i] < r[i][0] || p_[i] > r[i][1])
            throw exceptions::ParOutOfRange<Q,isDiscrete>
                {get_name(), i, p_[i], r[i], msgr};
}


template<typename T,
         bool isDiscrete>
auto Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_arithmetic_v<T>>>::
readin(
    const str& name,
    const std::filesystem::path& file,
    Msgr* msgr
)
{
    return Par<Q,isDiscrete> {name, file, msgr}();
}


template<typename T,
         bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_arithmetic_v<T>>>::
print(Msgr* msgr)
{
    if (msgr) {
        msgr->print_vector(get_name(), p_);
        return;
    }
    std::cout << "Messenger is not initialized for printing "
                << get_name() << std::endl;
}


template<typename T,
         bool isDiscrete>
std::vector<T> Par<std::vector<T>,
                   isDiscrete,
                   std::enable_if_t<std::is_arithmetic_v<T>>>::
operator()() const
{
    ASSERT(true, get_name());
    return p_;
}


template<typename T,
         bool isDiscrete>
T Par<std::vector<T>,
      isDiscrete,
      std::enable_if_t<std::is_arithmetic_v<T>>>::
operator[](const szt i) const
{
    ASSERT(isLoaded_, get_name());
    ASSERT(i<p_.size(), get_name());
    return p_[i];
}


template<typename T,
         bool isDiscrete>
void Par<std::vector<T>,
         isDiscrete,
         std::enable_if_t<std::is_arithmetic_v<T>>>::
initialize(str value)
{

    auto size_message = [&](const str& s, szt n) {
        return std::string
            {"Improper Config::" + get_name()+
             " initialization: " + s + " data: \nExpected are " +
             std::to_string(expectedSize_) + " elements, read " +
             std::to_string(n) + " elements"};
    };

    const str emp {" "};
    const str tab {"\t"};
    while (value.length()) {
        unsigned long e {value.find(emp)};
        if (e == str::npos) e = value.find(tab);
        if (e == str::npos) e = value.length();
        const str val {value.substr(0, e)};
        if (val.length() < 1)
            throw common::Exception
                {"Error in config file: Number of elelments in " + get_name() +
                 " is " + std::to_string(p_.size()) + " which is insufficient",
                 nullptr
                };
        T tmp;
        std::stringstream(val) >> tmp;
        p_.push_back(tmp);
        value.erase(0, e);
        value = common::trim(value);
        if (p_.size() > expectedSize_)
            throw common::Exception {size_message("Excessive", p_.size()), nullptr};
    }
    if (p_.size() != expectedSize_)
        throw common::Exception {size_message("Incorrect", p_.size()), nullptr};
}

}  // namespace utils::config::parameter

#endif // UTILS_CONFIG_PARAMETER_VECTORS_H
