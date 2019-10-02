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
 * \file scalars_strings.h
 * \brief Parameters of fundamental or string type.
 * \details Contains template partial specialization for classes encapsulating confuguration file
 * parameters of fundamental or string type.
 + \author Valerii Sukhorukov
 */

#ifndef UTILS_CONFIG_PARAMETER_SCALARS_STRINGS_H
#define UTILS_CONFIG_PARAMETER_SCALARS_STRINGS_H

#include <sstream>

#include "../../common/misc.h"
#include "../../common/msgr.h"
#include "base.h"
#include "../exceptions/scalars.h"

/// Library-wide.
namespace Utils {
/// Configuration module
namespace Config {
namespace Parameter {

using namespace Common;
using namespace Exceptions;

/**
* \brief Partial template specialization for for parameters of std fundamental scalar or std::string types.
* \tparam T Parameter type: must be std::is_fundamental.
* \tparam isDiscrete Specifies if the vector components accept discrete of continous values.
*/
template <typename T, bool isDiscrete>
class Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
												  std::is_same<T,std::string>::value>>
	: public Base<T> {

	using Q = T;
	using Base<T>::get_name;
	using Base<T>::isLoaded_;
	
	Q p_;		///< The parameter value.
	
public:	

	/**
	* \brief Constructor.
	* \param name Name of the parameter.
	*/
	explicit Par(const std::string& name);

	/**
	* \brief Constructor.
	* \param name Name of the parameter.
	* \param fname Name of the configuration file.
	* \param range Acceptable range of parameter values.
	* \param msgr \a Msgr used for the output.
	* \see Msgr
	*/
	explicit Par(
		 	const std::string& name,
		 	const std::string& fname,
		 	const std::vector<T>& range,
		 	Msgr* msgr=nullptr);

	/**
	* \brief Constructor.
	* \param name Name of the parameter.
	* \param fname Name of the configuration file.
	* \param range Acceptable range of parameter values.
	* \param msgr \a Msgr used for the output.
	* \see Msgr
	*/
	explicit Par(
		 	const std::string& name,
		 	const std::string& fname,
		 	const std::array<T,2>& range,
		 	Msgr* msgr=nullptr);

/*	void check_range( const std::array<T,2>& r, Msgr& msgr )
	{
		if (isDiscrete) {
			if (std::find(r.begin(), r.end(), p_) == r.end())
				throw Exceptions::ParOutOfRange<T,isDiscrete>(name(), p_, r, msgr);
		}
		else {
			XASSERT(r.size()==2, "size of r must be 2 for continuous parameters");
			if (p_<r[0] || p_>r[1])
				throw Exceptions::ParOutOfRange<T,isDiscrete>(name(), p_, r, msgr);
		}
	}
*/

	/**
	* \brief Check that the read in parameter value is within the range set by \p r.
	* \tparam W Acceptable range of parameter values.
	* \param r Acceptable range of parameter values.
	* \param msgr \a Msgr used for the output.
	* \see Msgr
	*/
	template <typename W>
	void check_range(const W& r,
					 Msgr* msgr=nullptr);

	/**
	* \brief Static function for reading a parameter without instantiating this class object.
	* \param name Name of the parameter.
	* \param fname Name of the configuration file.
	* \param range Acceptable range of parameter values.
	* \param msgr \a Msgr used for the output.
	* \return Parameter values (the whole vector).
	*/
	static auto readin(const std::string& name,
					   const std::string& fname,
					   const std::vector<Q>& range,
					   Msgr* msgr=nullptr);

	/**
	* \brief Print the the parameter to std::cout and logfile.
	* \param msgr \a Msgr used for the output.
	* \see Msgr
	*/
	void print(Msgr* msgr=nullptr) final;

	/**
	* \brief Set the field \a p_ to the value \p val.
	* \param val New value of the parameter.
	*/
	void set(const Q& val);

	/**
	* \brief The parameter value.
	* \return Parameter value.
	*/
	Q operator()() const;
	
private:

	/**
	* \brief Initialize the parameter from the config file.
	* \param value Value to search for.
	*/
	void initialize(std::string value) final;

};	

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
Par( const std::string& name )
	: Base<T> {name}
{}

template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
Par( const std::string& name,
	 const std::string& fname,
	 const std::vector<T>& range,
	 Msgr* msgr )
	: Base<T> {name}
{
	this->load(fname);
	check_range(range, msgr);
	print(msgr);
}

template <typename T, bool isDiscrete>
Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
Par( const std::string& name,
	 const std::string& fname,
	 const std::array<T,2>& range,
	 Msgr* msgr )
	: Base<T> {name}
{
	this->load(fname);
	check_range(range, msgr);
	print(msgr);
}

template <typename T, bool isDiscrete>
template <typename W>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
check_range( const W& r,
			 Msgr* msgr )
{
	if (!r.size()) return;		// use this case to omit string checkups

	if constexpr (isDiscrete) {
		if (std::find(r.begin(), r.end(), p_) == r.end())
			throw Exceptions::ParOutOfRange<T,isDiscrete>(get_name(), p_, r, msgr);
	}
	else {
		XASSERT(r.size()==2, "size of r must be 2 for continuous parameters");
		if (p_<r[0] || p_>r[1])
			throw Exceptions::ParOutOfRange<T,isDiscrete>(get_name(), p_, r, msgr);
	}
}

template <typename T, bool isDiscrete>
auto Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
readin( const std::string& name,
		const std::string& fname,
		const std::vector<Q>& range,
		Msgr* msgr )
{
	return Par<Q,isDiscrete> {name, fname, range, msgr}();
}

template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
print( Msgr* msgr )
{
	if (msgr)
		msgr->print(get_name()+" = "+STR(p_));
	else
		std::cout << get_name() << " = " << p_ << std::endl;
}

template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
set( const Q& val )
{
	isLoaded_ = true;
	p_ = val;
}

template <typename T, bool isDiscrete>
T Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
operator()() const
{
	XASSERT(isLoaded_, get_name());
	return p_;
}

template <typename T, bool isDiscrete>
void Par<T,isDiscrete, typename std::enable_if_t<std::is_fundamental<T>::value ||
											std::is_same<T,std::string>::value>>::
initialize( std::string value )
{
	std::stringstream(value) >> p_;
}

}	// namespace Parameter
}	// namespace Config
}   // namespace Utils

#endif // UTILS_CONFIG_PARAMETER_SCALARS_STRINGS_H
