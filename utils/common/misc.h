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


#ifndef UTILS_COMMON_MISC_H
#define UTILS_COMMON_MISC_H

#include <cmath>
#include <limits>
#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <numeric>
#include <fstream>
#include <memory>
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <sys/stat.h>
//#include <boost/filesystem.hpp>

#ifdef _DEBUG
	#define XASSERT(EX, msg) \
		(void)( (EX) || ( Utils::Common::assert_fun( #EX, __FILE__, __LINE__, msg ), 0 ) )
#else
	#define XASSERT(EX, msg)
#endif

//#include "../arrays/all.h"

namespace Utils {
namespace Common {

extern const std::string SLASH;

#define STR(x) std::to_string(x)

using ulong = unsigned long;
using uint = unsigned int;
using szt = std::size_t;

//template <typename T> using vec2 = std::vector<std::vector<T>>;
template <typename T, auto N> using arrvec = std::array<std::vector<T>,N>;
template <typename T, auto N> using vecarr = std::vector<std::array<T,N>>;
template <typename T> using vup = std::vector<std::unique_ptr<T>>;

template <typename T> constexpr T zero {static_cast<T>(0.L)};
template <typename T> constexpr T half {static_cast<T>(.5L)};
template <typename T> constexpr T thrd {static_cast<T>(1.L/3.L)};
template <typename T> constexpr T one {static_cast<T>(1.L)};
template <typename T> constexpr T two {static_cast<T>(2.L)};
template <typename T> constexpr T three {static_cast<T>(3.L)};
template <typename T> constexpr T four {static_cast<T>(4.L)};
template <typename T> constexpr T five {static_cast<T>(5.L)};
template <typename T> constexpr T six {static_cast<T>(6.L)};
template <typename T> constexpr T ten {static_cast<T>(10.L)};

template <typename T> constexpr T pi {static_cast<T>(3.1415926535897932384626433832795L)};
template <typename T> constexpr T twopi {two<T>*pi<T>};
template <typename T> constexpr T halfpi {half<T>*pi<T>};
template <typename T> constexpr T sqrtPI {static_cast<T>(1.7724538509055160272981674833411L)};

template <typename T, typename Enabler = void> constexpr T sqrt2PI;
template <typename T> constexpr T sqrt2PI<T,std::enable_if_t<std::is_arithmetic<T>::value>> {std::pow(twopi<T>, half<T>)};

template <typename T, typename Enabler = void> constexpr T EPS;
template <typename T> constexpr const T EPS<T,std::enable_if_t<std::is_fundamental<T>::value>> {std::numeric_limits<T>::epsilon()};

template <typename T, typename Enabler = void> constexpr T huge;
template <typename T> constexpr const T huge<T,std::enable_if_t<std::is_fundamental<T>::value>>
	{std::numeric_limits<T>::has_infinity
   ? std::numeric_limits<T>::infinity()
   : std::numeric_limits<T>::max()};
template <typename T, typename Enabler = void> constexpr T MAX;
template <typename T> constexpr T MAX<T,std::enable_if_t<std::is_fundamental<T>::value>> = std::numeric_limits<T>::max();
template <typename T> constexpr T MAX<T,std::enable_if_t<std::is_same<T,std::string>::value>> {""};
template <typename T, typename Enabler = void> constexpr T MIN;
template <typename T> constexpr T MIN<T,std::enable_if_t<std::is_fundamental<T>::value>> {std::numeric_limits<T>::min()};
template <typename T> constexpr T MIN<T,std::enable_if_t<std::is_same<T,std::string>::value>> {""};

template <typename T> constexpr T INF {std::numeric_limits<T>::infinity()};

template <typename T, auto N> constexpr
auto filled_array( const T val )
{
	std::array<T,N> a {};
	for (auto& o : a)
		o = val;
//		a.fill(val);	// needs c++20 to be constexpr
	return a;
}

template <auto N> constexpr std::array<bool,N> falses {filled_array<bool,N>(false)};
template <auto N> constexpr std::array<bool,N> trues {filled_array<bool,N>(true)};
template <typename T, auto N> constexpr std::array<T,N> zeros {filled_array<T,N>(zero<T>)};
template <typename T, auto N> constexpr std::array<T,N> ones {filled_array<T,N>(one<T>)};
template <typename T, auto N> constexpr std::array<T,N> hundreds {filled_array<T,N>(static_cast<T>(100.L))};
template <typename T, auto N> constexpr std::array<T,N> huges {filled_array<T,N>(huge<T>)};
template <typename T, auto N> constexpr std::array<T,N> mhuges {filled_array<T,N>(-huge<T>)};
constexpr std::array<bool,2> bools {{false, true}};
template <typename T> constexpr std::array<T,2> zeroone {{zero<T>, one<T>}};
template <typename T> constexpr std::array<T,2> zerohuge {{zero<T>, huge<T>}};
template <typename T> constexpr std::array<T,2> onehuge {{one<T>, huge<T>}};
template <typename T> constexpr std::array<T,2> moneone {{-one<T>, one<T>}};
template <typename T> constexpr std::array<T,2> mhugehuge {{-huge<T>, huge<T>}};
template <auto N> constexpr const std::vector rangeBools {falses<N>, trues<N>};
template <typename T, auto N> constexpr auto vecarr0H() { return vecarr<T,N>{zeros<T,N>, huges<T,N>}; }
template <typename T, auto N> constexpr auto vecarr01() { return vecarr<T,N>{zeros<T,N>, ones<T,N>}; }
template <szt N> 			 constexpr auto vecarrFT() { return vecarr<bool,N>{falses<N>, trues<N>}; }
//template <typename T, auto N> constexpr const std::vector<T,N> zeroshuges {zeros<real,N>, huges<real,N>};

std::string operator"" _str (long double number);
std::string operator"" _str (unsigned long long number);

template <typename T, typename Q, T (Q::* P)() const> // member function pointer parameter
struct Adder {
    int operator()(const T& i, const Q& o) const {
        return (o.*P)() + i;
    }
};
/*
template <typename T>
class vec2 {
	std::vector<std::vector<T>> v;
public:
	vec2() = default;
	vec2( const szt x, const szt y, const T ini=zero<T> )
	{
		v.resize(x);
		for (auto& vv : v)
			vv.resize(y, ini);
	}
	constexpr std::vector<T> operator[]( szt i ) const noexcept {
		return v[i];
	}
	std::vector<T>& operator[]( szt i ) noexcept {
		return v[i];
	}
	auto size() const noexcept { return v.size(); }
	void resize(const szt n) noexcept { v.resize(n); }
	auto begin() noexcept { return v.begin(); }
	auto end() noexcept {return v.end(); }
	auto clear() noexcept { v.clear(); }
	auto push_back(const std::vector<T>& u) noexcept { v.push_back(u); }
	void fill( const T val )
	{
		for (auto& o : v)
			for (auto& oo : o)
				oo = val;
	}

};*/
long long assert_fun(
	const char* EX,
	const char *file,
	int line,
	const std::string& msg );

template <typename T> using vec2 = std::vector<std::vector<T>>;
template <typename T> using vec3 = std::vector<vec2<T>>;
template <typename T> using vec4 = std::vector<vec3<T>>;

namespace Vec2 {
template <typename T1, typename T2> inline
vec2<T1> array_like( const vec2<T2>& as )
{
	vec2<T1> me(as.size());
	for (szt i=0; i<as.size(); i++)
		me[i].resize(as[i].size());
	return me;
}

template <typename T> inline
vec2<T> make( const szt x, const szt y, const T ini=zero<T> )
{
	vec2<T> v(x);
	for (auto& vv : v)
		vv.resize(y, ini);
	return v;
}

template <typename T> inline
szt size( const vec2<T>& v )
{
	szt s {};
	for (auto& vv : v)
		s += vv.size();
	return s;
}

template <typename T>
void add_scalar( const T d, vec2<T>& p )
{
	for (auto& o : p)
		for (auto& oo : o)
			oo += d;
}
template <typename T> inline
void fill( vec2<T>& v, const T val )
{
	for (auto& o : v)
		for (auto& oo : o)
			oo = val;
}

}  // namespace Vec2

namespace Vec3 {

template <typename T> inline
vec3<T> make( const szt x, const szt y, const szt z, const T ini=zero<T> )
{
	vec3<T> v(x);
	for (auto& vv : v) {
		vv.resize(y);
		for (auto& vvv : vv)
			vvv.resize(z, ini);
	}
	return v;
}

}	// namespace Vec3

template <typename T> inline
void partial_sum( T const* u, T* v, const szt from, const szt num ) noexcept
{
	const auto f {static_cast<int>(from)};
	const auto n {static_cast<int>(num)};

	v[f] = u[f];
	if (f) v[f] += v[f-1];

	for (int i=f+1; i<f+n; i++)
		v[i] = v[i-1] + u[i];
}

template <typename T> inline
T avg(const std::vector<T>& v)
{
    return std::accumulate(v.begin(), v.end(), zero<T>) / v.size();
}
double avg(std::vector<int> const& v);

template <typename T> inline
T var( const std::vector<T>& n )
{
	T v {};
	T mn {avg(n)};
	for (const auto& o : n)
		v += (o - mn) * (o - mn);
	return v / n.size();
}

std::string trim(
	const std::string& str,
	const std::string& whitespace = " "
);

template <typename T> constexpr
std::vector<T> exp_num( const T b, const T r, const T dx ) noexcept
{
	uint n {r / dx};
	std::vector<T> x(n), q(n);

	x[0] = zero<T>;
	for (uint i=1; i<n; i++)
		x[i] = x[i-1] + dx;

	const auto c {b / (std::exp( b * r) - one<T>)};

	for (uint i=0; i<n; i++)
		q[i] = c * std::exp(b * x[i]);

	return q;
}

template <typename T> constexpr
T sigmoid_decay( const T x, const T x0, const T steepness ) noexcept
{
	return one<T> / (one<T> + std::exp((x - x0) * steepness));
}

// returns how many elements in b /= 0 putting their indices to j
template <typename T> 
szt find( const std::vector<T>& b, 
		  std::vector<szt>& j )	noexcept
{
	j.clear();
	for (szt i=0; i<b.size(); i++) 
		if (b[i] != zero<T>)
			j.push_back(i);
	return j.size();
}

template <typename T>
bool find_firstFromFront( const std::vector<T>& b, const T w, szt& j )
{
	if (!b.size())
		return false;

	for (szt i=0; i<b.size(); i++)
		if (b[i] == w) {
			j = i;
			return true;
		}
	return false;
}

template <typename T>
bool find_firstFromBack( const std::vector<T>& b, const T w, szt& j )
{
	if (!b.size())
		return false;

	for (long i=b.size()-1; i>=0; i--)
		if (b[i] == w) {
			j = static_cast<szt>(i);
			return true;
		}
	return false;
}

template <typename T> constexpr
szt index_min( const std::vector<T>& v ) noexcept
{
	szt k {};
	auto m {v[0]};
	for (szt i=1; i<v.size(); i++)
		if (v[i] < m) {
			k = i;
			m = v[i];
		}
	return k;
}

// standard Gaussian function
template <typename T> constexpr
T gaussian( const T x ) noexcept
{
	return std::exp(-x*x/two<T>)/std::sqrt(twopi<T>);
}

template <typename T> constexpr
T gaussian( const T x, const T var ) noexcept
{
	return std::exp(-x*x/(two<T>*var))/std::sqrt(twopi<T>*var);
}

template <typename T> constexpr
T gaussian( const T x, const T mean, const T var ) noexcept
{
	return std::exp(-(x-mean)*(x-mean)/(two<T>*var))/std::sqrt(twopi<T>*var);
}

template <typename T>
T gaussian_fun( T x, T mean, T sigma )
{
	return one<T> / std::sqrt(twopi<T>) * std::exp(-(x-mean)*(x-mean)/(sigma*sigma)/two<T>);
}

template <typename T> constexpr
auto gr2rad( T grad )
{
	return grad*pi<T>/static_cast<T>(180.L);
}

template <typename T> constexpr
auto rad2gr( T rad )
{
	return rad*static_cast<T>(180.L)/pi<T>;
}

template <auto K>
std::string pad_zeros( const szt n )
{
	static_assert(K > 1 && K < 7, "Padding is only supprted for lengths between 2 and 6 inclusive");

	if constexpr (K == 2)
		return n<10 ? "0" : "" + STR(n);
	else if constexpr (K == 3)
		return n<100 ? n<10 ? "00" : "0" : "" + STR(n);
	else if constexpr (K == 4)
		return n<1000 ? n<100 ? n<10 ? "000" : "00" : "0" : "" + STR(n);
	else if constexpr (K == 5)
		return n<10000 ? n<1000 ? n<100 ? n<10 ? "0000" : "000" : "00" : "0" : "" + STR(n);
	else if constexpr (K == 6)
		return n<100000 ? n<10000 ? n<1000 ? n<100 ? n<10 ? "00000" : "0000" : "000" : "00" : "0" : "" + STR(n);
};

bool file_exists( const std::string& name );
//bool fileExists( const std::string& name );
bool directory_exists( const std::string& pathstrconst );
//void check_directory(const std::string& s );
void copy_text_file(const std::string& fname1, const std::string& fname2);

struct StopWatch {

	struct Instance {

		std::chrono::time_point<std::chrono::system_clock> h;
		std::time_t c;
		std::string str;

		void operator()() {
			h = std::chrono::system_clock::now();
			c = std::chrono::system_clock::to_time_t(h);
			str = std::string(ctime(&c));
		}
	};

	Instance start;
	Instance stop;

	std::string duration() {
		diff = stop.h - start.h;
		return STR(diff.count());
	}

private:

	std::chrono::duration<double> diff;

};

template <typename F, typename S>
std::vector<F> map_first_as_vec(const std::map<F,S>& ct )
{
	std::vector<F> v;
	for (const auto& o : ct)
		v.emplace_back(o.first);
	return v;
}

template <typename F, typename S>
std::vector<S> map_second_as_vec(const std::map<F,S>& ct )
{
	std::vector<S> v;
	for (const auto& o : ct)
		v.emplace_back(o.second);
	return v;
}

template <typename T>
void remove_vector_element( std::vector<T>& v, const T& a )
{
	for (auto i=v.begin(); i!=v.end(); i++)
		if (a[0] == (*i)[0] &&
			a[1] == (*i)[1] ) {
			v.erase(i);
			return;
		}
	XASSERT(false, "Error in remove_vector_element: element not found");
}


#define ANSI_RESET       "\x1B[0m"
#define ANSI_FG_BLACK    "\x1b[30m"
#define ANSI_FG_RED      "\x1B[31m"
#define ANSI_FG_GREEN    "\x1B[32m"
#define ANSI_FG_YELLOW   "\x1B[33m"
#define ANSI_FG_BLUE     "\x1B[34m"
#define ANSI_FG_MAGENTA  "\x1B[35m"
#define ANSI_FG_CYAN     "\x1B[36m"
#define ANSI_FG_WHITE    "\x1B[37m"
#define ANSI_BG_RED      "\x1b[41m"
#define ANSI_BG_GREEN    "\x1b[42m"
#define ANSI_BG_YELLOW   "\x1b[43m"
#define ANSI_BG_BLUE     "\x1b[44m"
#define ANSI_BG_MAGENTA  "\x1b[45m"
#define ANSI_BG_CYAN     "\x1b[46m"
#define ANSI_BG_WHITE    "\x1b[47m"
#define ANSI_BOLD_ON     "\x1b[1m"			
#define ANSI_BOLD_OFF    "\x1b[22m"		
#define ANSI_INVERSE_ON  "\x1b[7m"


} 	// namespace Common
}	// namespace Utils


#endif // UTILS_COMMON_MISC_H
