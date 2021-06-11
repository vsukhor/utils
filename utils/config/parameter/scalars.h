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
 * \file scalars.h
 * \brief Parameters of arithmetic type.
 * \details Contains template partial specialization for classes encapsulating
 * confuguration file parameters of arithmetic type.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_SCALARS_H
#define UTILS_CONFIG_PARAMETER_SCALARS_H

#include <filesystem>
#include <sstream>
#include <string>

#include "../../common/misc.h"
#include "../../common/msgr.h"
#include "../exceptions/scalars.h"
#include "base.h"

/// Parameter namespace.
namespace utils::config::parameter {

/**
 * Parameters of arithmetic scalar types.
 * Partial template specialization for for parameters of arithmetic
 * scalar types.
 * \tparam T Parameter type: must be std::is_arithmetic continous values.
 */
template <typename T, bool isDiscrete>
class Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>
    : public Base<T> {

    using Q = T;
    using Base<T>::check_name;
    using Base<T>::isLoaded_;
    
    Q p_;  ///< The parameter value.
    
public:    

    using Base<T>::get_name;

    /**
     * Constructor.
     * \param name Name of the parameter.
     */
    explicit Par(const std::string& name);

    /**
     * Constructor.
     * \param name Name of the parameter.
     * \param file Configuration file.
     * \param range Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    explicit Par(
        const std::string& name,
        const std::filesystem::path& file,
        const std::vector<T>& range,
        Msgr* msgr=nullptr
    );

    /**
     * Constructor.
     * \param name Name of the parameter.
     * \param file Configuration file.
     * \param range Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    explicit Par(
        const std::string& name,
        const std::filesystem::path& file,
        const std::array<T,2>& range,
        Msgr* msgr=nullptr
    );

    /**
     * Check that the read in parameter value is within the range set by \p r.
     * \tparam W Acceptable range of parameter values.
     * \param r Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    template <typename W>
    void check_range(const W& r,
                     Msgr* msgr=nullptr);

    /**
     * Read a parameter without instantiating.
     * Static function for reading a parameter without instantiating
     * this class object.
     * \param name Name of the parameter.
     * \param file Configuration file.
     * \param range Acceptable range of parameter values.
     * \param msgr \a Msgr used for the output.
     * \return Parameter values (the whole vector).
     */
    static auto readin(
        const std::string& name,
        const std::filesystem::path& file,
        const std::vector<Q>& range,
        Msgr* msgr=nullptr
    );

    /**
     * Print the the parameter to std::cout and logfile.
     * \param msgr \a Msgr used for the output.
     * \see Msgr
     */
    void print(Msgr* msgr=nullptr) final;

    /**
     * Set the field \a p_ to the value \p val.
     * \param val New value of the parameter.
     */
    void set(const Q& val);

    /**
     * The parameter value.
     * \return Parameter value.
     */
    Q operator()() const;
    
private:

    /**
     * Initialize the parameter from the config file.
     * \param value Value to search for.
     */
    void initialize(std::string value) final;
};    


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
Par( const std::string& name )
    : Base<T> {check_name(name)}
{}


template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
Par( const std::string& name,
     const std::filesystem::path& file,
     const std::vector<T>& range,
     Msgr* msgr )
    : Base<T> {check_name(name)}
{
    this->load(file);
    try {
        check_range(range, msgr);
    }
    catch (const exceptions::ParOutOfRange<T,isDiscrete>&) {
        std::exit(EXIT_FAILURE);
    }
    print(msgr);
}


template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
Par( const std::string& name,
     const std::filesystem::path& file,
     const std::array<T,2>& range,
     Msgr* msgr )
    : Base<T> {check_name(name)}
{
    this->load(file);
    try {
        check_range(range, msgr);
    }
    catch (const exceptions::ParOutOfRange<T,isDiscrete>&) {
        std::exit(EXIT_FAILURE);
    }
    print(msgr);
}


template <typename T, bool isDiscrete>
template <typename W>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
check_range( const W& r,
             Msgr* msgr )
{
    if (!r.size()) return;        // use this case to omit string checkups

    if constexpr (isDiscrete) {
        if (std::find(r.begin(), r.end(), p_) == r.end())
            throw exceptions::ParOutOfRange<T,isDiscrete>(get_name(), p_, r, msgr);
    }
    else {
        XASSERT(r.size() == 2, "size of r must be 2 for continuous parameters");
        if (p_ < r[0] || p_ > r[1])
            throw exceptions::ParOutOfRange<T,isDiscrete>(get_name(), p_, r, msgr);
    }
}


template <typename T, bool isDiscrete>
auto Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
readin( const std::string& name,
        const std::filesystem::path& file,
        const std::vector<Q>& range,
        Msgr* msgr )
{
    return Par<Q,isDiscrete> {name, file, range, msgr}();
}


template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
print( Msgr* msgr )
{
    if (msgr != nullptr)
        msgr->print(get_name(), " = ", p_);
    else
        std::cout << get_name() << " = " << p_ << std::endl;
}


template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
set( const Q& val )
{
    isLoaded_ = true;
    p_ = val;
}


template <typename T, bool isDiscrete>
T Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
operator()() const
{
    XASSERT(isLoaded_, get_name());
    return p_;
}


template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_arithmetic_v<T>>>::
initialize( std::string value )
{
    std::stringstream(value) >> p_;
}

}  // namespace utils::config::parameter

#endif // UTILS_CONFIG_PARAMETER_SCALARS_H
