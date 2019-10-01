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
* \file array4.h
* \brief Four-element arrays.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_ARRAYS_ARRAY4_H
#define UTILS_ARRAYS_ARRAY4_H

#include <fstream>

#include "../common/misc.h"

/// Library-wide.
namespace Utils {
/// Custom arrays.
namespace Arrays {

/// \brief Four-element arrays.
/// \details This class specializes array template for four-element array of arithmetic types.
/// Implements convenient arithmetics as well as some functionaity
/// commonly used in 4-dimensional geometric applications.
/// \tparam T Type of the elements.
template <typename T>
class array<4,T,std::enable_if_t<std::is_arithmetic<T>::value>> {

T n[4] = {};
 
public:
array( T m=static_cast<T>(0) ) noexcept
	: n {m, m, m, m}
{}
array( T n0, T n1, T n2, T n3 ) noexcept
	: n {n0, n1, n2, n3}
{}
array( const array<2,T>& n1, const array<2,T>& n2 ) noexcept
	: n {n1[0], n1[1], n2[0], n2[1]}
{}
array( T n1, const array<3,T>& n2 ) noexcept
	: n {n1, n2[0], n2[1], n2[2]}
{}
array( const array<3,T>& n1, T n2 ) noexcept
	: n {n1[0], n1[1], n1[2], n2}
{}
array( const array& p ) noexcept
	: n {p[0], p[1], p[2], p[3]}
{}
array( const std::array<T,4>& p ) noexcept
	: n {p[0], p[1], p[2], p[3]}
{}
template<typename Q>
constexpr array<4,Q> cast_static() const noexcept {
	return {static_cast<Q>(n[0]),
			static_cast<Q>(n[1]),
			static_cast<Q>(n[2]),
			static_cast<Q>(n[3])};
}
constexpr array operator=( const array& p ) noexcept {
	n[0] = p[0];
	n[1] = p[1];
	n[2] = p[2];
	n[3] = p[3];
	return *this;
}
constexpr array operator=( const std::array<T,4>& p ) noexcept {
	n[0] = p[0];
	n[1] = p[1];
	n[2] = p[2];
	n[3] = p[3];
	return *this;
}
constexpr array operator=( const T p[] ) noexcept {
	n[0] = p[0];
	n[1] = p[1];
	n[2] = p[2];
	n[3] = p[3];
	return *this;
}
constexpr array operator=( T p ) noexcept {
	n[0] = p;
	n[1] = p;
	n[2] = p;
	n[3] = p;
	return *this;
}
constexpr array<2,T> operator()( int i1, int i2 ) const noexcept {
	XASSERT(i1 >= 0, "");
	XASSERT(i1 <  4, "");
	XASSERT(i2 >= 0, "");
	XASSERT(i2 <  4, "");
	return {n[i1], n[i2]};
}
constexpr array<3,T> operator()( int i1, int i2, int i3 ) const noexcept {
	XASSERT(i1 >= 0, "");
	XASSERT(i1 <  4, "");
	XASSERT(i2 >= 0, "");
	XASSERT(i2 <  4, "");
	XASSERT(i3 >= 0, "");
	XASSERT(i3 <  4, "");
	return {n[i1], n[i2], n[i3]};
}
constexpr array operator+( const array& p ) const noexcept {
	return { n[0] + p[0],
			 n[1] + p[1],
			 n[2] + p[2],
			 n[3] + p[3] };
}
constexpr array operator+( const T p[] ) const noexcept {
	return { n[0] + p[0],
			 n[1] + p[1],
			 n[2] + p[2],
			 n[3] + p[3] };
}
constexpr array& operator+=( const array& p ) noexcept {
	n[0] += p[0];
	n[1] += p[1];
	n[2] += p[2];
	n[3] += p[3];
	return *this;
}
constexpr array& operator+=( const T p[] ) noexcept {
	n[0] += p[0];
	n[1] += p[1];
	n[2] += p[2];
	n[3] += p[3];
	return *this;
}
constexpr array operator+( T p ) const noexcept {
	return { n[0] + p,
			 n[1] + p,
			 n[2] + p,
			 n[3] + p };
}
constexpr array operator-() const noexcept {
	array q;
	q[0] = -n[0];
	q[1] = -n[1];
	q[2] = -n[2];
	q[3] = -n[3];
	return q;
}
constexpr array operator-( const array& p ) const noexcept {
	return { n[0] - p[0],
			 n[1] - p[1],
			 n[2] - p[2],
			 n[3] - p[3] };
}
constexpr array operator-( const T p[] ) const noexcept {
	return { n[0] - p[0],
			 n[1] - p[1],
			 n[2] - p[2],
			 n[3] - p[3] };
}
constexpr array& operator-=( const array& p ) noexcept {
	n[0] -= p[0];
	n[1] -= p[1];
	n[2] -= p[2];
	n[3] -= p[3];
	return *this;
}
constexpr array& operator-=( const T p[] ) noexcept {
	n[0] -= p[0];
	n[1] -= p[1];
	n[2] -= p[2];
	n[3] -= p[3];
	return *this;
}
constexpr array operator-( T p ) const noexcept {
	return { n[0] - p,
			 n[1] - p,
			 n[2] - p,
			 n[3] - p };
}
constexpr array operator*( const array& p ) const noexcept {
	return { n[0] * p[0],
			 n[1] * p[1],
			 n[2] * p[2],
			 n[3] * p[3] };
}
constexpr array operator*( const T p[] ) const noexcept {
	return { n[0] * p[0],
			 n[1] * p[1],
			 n[2] * p[2],
			 n[3] * p[3] };
}
constexpr array& operator*=( const array& p ) noexcept {
	n[0] *= p[0];
	n[1] *= p[1];
	n[2] *= p[2];
	n[3] *= p[3];
	return *this;
}
constexpr array& operator*=( const T p[] ) noexcept {
	n[0] *= p[0];
	n[1] *= p[1];
	n[2] *= p[2];
	n[3] *= p[3];
	return *this;
}
constexpr array operator*( T p ) const noexcept {
	return { n[0] * p,
			 n[1] * p,
			 n[2] * p,
			 n[3] * p };
}
constexpr array operator/( const array& p ) const noexcept {
	return { n[0] / p[0],
			 n[1] / p[1],
			 n[2] / p[2],
			 n[3] / p[3] };
}
constexpr array operator/( const T p[] ) const noexcept {
	return { n[0] / p[0],
			 n[1] / p[1],
			 n[2] / p[2],
			 n[3] / p[3] };
}
constexpr array& operator/=( const array& p ) noexcept {
	n[0] /= p[0];
	n[1] /= p[1];
	n[2] /= p[2];
	n[3] /= p[3];
	return *this;
}
constexpr array& operator/=( const T p[] ) noexcept {
	n[0] /= p[0];
	n[1] /= p[1];
	n[2] /= p[2];
	n[3] /= p[3];
	return *this;
}
constexpr array operator/( T p ) const noexcept {
	return { n[0] / p,
			 n[1] / p,
			 n[2] / p,
			 n[3] / p };
}
constexpr bool operator==( const array& p ) const noexcept {
	return n[0] == p[0] && 
		   n[1] == p[1] && 
		   n[2] == p[2] && 
		   n[3] == p[3];
}
constexpr bool operator==( const std::array<T,4>& p ) const noexcept {
	return n[0] == p[0] && 
		   n[1] == p[1] && 
		   n[2] == p[2] && 
		   n[3] == p[3];
}
constexpr bool operator==( const T p[] ) const noexcept {
	return n[0] == p[0] && 
		   n[1] == p[1] && 
		   n[2] == p[2] && 
		   n[3] == p[3];
}
constexpr bool operator==( T p ) const noexcept {
	return n[0] == p && 
		   n[1] == p && 
		   n[2] == p && 
		   n[3] == p;
}
constexpr bool operator!=( const array& p ) const noexcept {
	return n[0] != p[0] || 
		   n[1] != p[1] || 
		   n[2] != p[2] || 
		   n[3] != p[3];
}
constexpr bool operator!=( const T p[] ) const noexcept {
	return n[0] != p[0] || 
		   n[1] != p[1] || 
		   n[2] != p[2] || 
		   n[3] != p[3];
}
constexpr bool operator!=( T p ) const noexcept {
	return n[0] != p || 
		   n[1] != p || 
		   n[2] != p || 
		   n[3] != p;
}
constexpr bool operator<( const array& p ) const noexcept {
	return n[0] < p[0] && 
		   n[1] < p[1] && 
		   n[2] < p[2] && 
		   n[3] < p[3];
}
constexpr bool operator<( T p ) const noexcept {
	return n[0] < p && 
		   n[1] < p && 
		   n[2] < p && 
		   n[3] < p;
}
constexpr bool operator<=( const array& p ) const noexcept {
	return n[0] <= p[0] && 
		   n[1] <= p[1] && 
		   n[2] <= p[2] && 
		   n[3] <= p[3];
}
constexpr bool operator<=( T p ) const noexcept {
	return n[0] <= p && 
		   n[1] <= p && 
		   n[2] <= p && 
		   n[3] <= p;
}
constexpr bool operator>( const array& p ) const noexcept {
	return n[0] > p[0] && 
		   n[1] > p[1] && 
		   n[2] > p[2] && 
		   n[3] > p[3];
}
constexpr bool operator>( T p ) const noexcept {
	return n[0] > p && 
		   n[1] > p && 
		   n[2] > p && 
		   n[3] > p;
}
constexpr bool operator>=( const array& p ) const noexcept {
	return n[0] >= p[0] && 
		   n[1] >= p[1] && 
		   n[2] >= p[2] && 
		   n[3] >= p[3];
}
constexpr bool operator>=( T p ) const noexcept {
	return n[0] >= p && 
		   n[1] >= p && 
		   n[2] >= p && 
		   n[3] >= p;
}
constexpr T operator[]( int i ) const noexcept {
	return n[i];
}
constexpr T& operator[]( int i ) noexcept {
	return n[i];
}
constexpr bool contains( T p ) const noexcept {
	return n[0] == p || 
		   n[1] == p || 
		   n[2] == p || 
		   n[3] == p;
}
constexpr void reflect() noexcept {
	T temp = n[0];
	n[0] = n[3];
	n[3] = temp;
	temp = n[1];
	n[1] = n[2];
	n[2] = temp;
}
constexpr T dotpr() const noexcept {
	return ( n[0]*n[0] + 
			 n[1]*n[1] + 
			 n[2]*n[2] + 
			 n[3]*n[3] );
}
constexpr T dotpr( const array& a ) const noexcept {
	return ( n[0]*a.n[0] + 
			 n[1]*a.n[1] + 
			 n[2]*a.n[2] + 
			 n[3]*a.n[3] );
}
static constexpr T dotpr( const array& a1,
						  const array& a2 ) noexcept {
	return ( a1.n[0]*a2.n[0] + 
			 a1.n[1]*a2.n[1] + 
			 a1.n[2]*a2.n[2] + 
			 a1.n[3]*a2.n[3] );
}
constexpr T norm() const noexcept {
	return std::sqrt(dotpr());
}
constexpr array unitv() const noexcept {
	return *this / norm();
}
constexpr T scaProjection( const array& b ) const noexcept {			// scalar projection of *this onto array b
	return dotpr(b) / b.norm();
}
constexpr array vecProjection( const array& b ) const noexcept {		// std::vector projection of *this onto array b
	return b.unitv() * scaProjection(b);
}
constexpr int find( T p ) noexcept {
	return ( p == n[0] ) ? 0 : 
		 ( ( p == n[1] ) ? 1 : 
		 ( ( p == n[2] ) ? 2 : 
		 ( ( p == n[3] ) ? 3 : -1 ) ) );
}
array<2,int> find( const array<2,T>& p ) noexcept {
	return ( p[0] == n[0] && p[1] == n[0] ) ? array<2,int>(0, 1) :
		 ( ( p == n[1] ) ? 1 : 
		 ( ( p == n[2] ) ? 2 : 
		 ( ( p == n[3] ) ? 3 : -1 ) ) );
}
array<3,T> other_than( T p ) noexcept {
	return ( p == n[0] ) ? array<3,T>( n[1], n[2], n[3] ) :
		 ( ( p == n[1] ) ? array<3,T>( n[0], n[2], n[3] ) :
		 ( ( p == n[2] ) ? array<3,T>( n[0], n[1], n[3] ) :
		 ( ( p == n[3] ) ? array<3,T>( n[0], n[1], n[2] ) : array<3,T>( -1 ) ) ) );
}	
array<2,T> other_than( const array<2,T>& p ) noexcept {
	return ( p == array<2,T>( n[0], n[1] ) || p == array<2,T>( n[1], n[0] ) ) ? array<2,T>( n[2], n[3] ) :
		 ( ( p == array<2,T>( n[0], n[2] ) || p == array<2,T>( n[2], n[0] ) ) ? array<2,T>( n[1], n[3] ) :
		 ( ( p == array<2,T>( n[0], n[3] ) || p == array<2,T>( n[3], n[0] ) ) ? array<2,T>( n[1], n[2] ) :
		 ( ( p == array<2,T>( n[1], n[2] ) || p == array<2,T>( n[2], n[1] ) ) ? array<2,T>( n[0], n[3] ) :
		 ( ( p == array<2,T>( n[1], n[3] ) || p == array<2,T>( n[3], n[1] ) ) ? array<2,T>( n[0], n[2] ) :
		 ( ( p == array<2,T>( n[2], n[3] ) || p == array<2,T>( n[3], n[2] ) ) ? array<2,T>( n[0], n[1] ) : array<2,T>( -1 ) ) ) ) ) );
}	
T other_than( const array<3,T>& p ) noexcept {
	return ( p == array<3,T>( n[0], n[1], n[2] ) || p == array<3,T>( n[0], n[2], n[1] ) || p == array<3,T>( n[1], n[0], n[2] ) || p == array<3,T>( n[1], n[2], n[0] ) || p == array<3,T>( n[2], n[0], n[1] ) || p == array<3,T>( n[2], n[1], n[0] ) ) ? n[3] :
		 ( ( p == array<3,T>( n[0], n[1], n[3] ) || p == array<3,T>( n[0], n[3], n[1] ) || p == array<3,T>( n[1], n[0], n[3] ) || p == array<3,T>( n[1], n[3], n[0] ) || p == array<3,T>( n[3], n[0], n[1] ) || p == array<3,T>( n[3], n[1], n[0] ) ) ? n[2] :
		 ( ( p == array<3,T>( n[0], n[3], n[2] ) || p == array<3,T>( n[0], n[2], n[3] ) || p == array<3,T>( n[3], n[0], n[2] ) || p == array<3,T>( n[3], n[2], n[0] ) || p == array<3,T>( n[2], n[0], n[3] ) || p == array<3,T>( n[2], n[3], n[0] ) ) ? n[1] :
		 ( ( p == array<3,T>( n[3], n[1], n[2] ) || p == array<3,T>( n[3], n[2], n[1] ) || p == array<3,T>( n[1], n[3], n[2] ) || p == array<3,T>( n[1], n[2], n[3] ) || p == array<3,T>( n[2], n[3], n[1] ) || p == array<3,T>( n[2], n[1], n[3] ) ) ? n[0] : -1 ) ) );
}	
constexpr T sum() const noexcept {
	return n[0] + n[1] + n[2] + n[3];
}
void print( std::ostream& os, bool end ) const noexcept {
	os << n[0] << " " << n[1] << " " << n[2] << " " << n[3];
	if (end) os << std::endl;
}
static constexpr array sum( const array a[],
				   const size_t& i1, 
				   const size_t& i2 ) noexcept {
	array res{};
	for (size_t i=i1; i<=i2; i++)
		res += a[i];
	return res;
}
constexpr static array mean( const std::vector<array>& a ) noexcept {
	return sum(a) / a.size();
}
static constexpr array mean( const array a[], size_t i1, size_t i2 ) noexcept {
	return sum(a, i1, i2) / (i2 - i1 + 1);
}
void read( std::ifstream& ist ) noexcept {
	ist.read( reinterpret_cast<char*>(&n[0]), sizeof(T));
	ist.read( reinterpret_cast<char*>(&n[1]), sizeof(T));
	ist.read( reinterpret_cast<char*>(&n[2]), sizeof(T));
	ist.read( reinterpret_cast<char*>(&n[3]), sizeof(T));
}
void write( std::ofstream& ost ) const noexcept {
	ost.write( reinterpret_cast<const char*>(&n[0]), sizeof(T));
	ost.write( reinterpret_cast<const char*>(&n[1]), sizeof(T));
	ost.write( reinterpret_cast<const char*>(&n[2]), sizeof(T));
	ost.write( reinterpret_cast<const char*>(&n[3]), sizeof(T));
}
};

}	// namespace Arrays
}	// namespace Utils

#endif // UTILS_ARRAYS_ARRAY4_H
