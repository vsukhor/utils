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
#ifndef UTILS_COMMON_GEOMETRIC_FUNCTIONS
#define UTILS_COMMON_GEOMETRIC_FUNCTIONS

#include "misc.h"
#include "msgr.h"
#include "../arrays/all.h"

namespace Utils {
namespace Common {

using namespace Arrays;

/**
* \class Geometric geometric_functions.h
* \brief A realtively loose collection of geometry-related static functions.
* \tparam Floating point type.
*/
template <typename T>
class Geometric {

public:

	/// Make sure that the template parameter is a floating type.
	static_assert(std::is_floating_point<T>::value,
				  "Class Geometric can only be instantiated with floating point types");


	// Elliptic shapes +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// Get a point on an ellipse centered at zero.
	static A2<T> ellipse(const T alpha,		///< Angular coordinate.
						 const A2<T>& ab	///< Dimensions of ellipse semi-major axes.
						 ) noexcept;

	/// Find out if point \p p is inside an ellipsoid given by dimensions \p e of the semi-major axes.
	static constexpr bool is_inside_ellipsoid(const A3<T>& p,
											  const A3<T>& e
											 ) noexcept;
	
	/// Calculate area of an ellipse.
	/// Its semi-major axes are \p a and \p b.
	static constexpr T ellipse_area(const T a,
									const T b
									) noexcept;

	/// \brief Calculate volume of an ellipsoid.
	/// Its semi-major axes are \p a, \p b and \p c..
	static constexpr T ellipsoid_vol(const T a,
									 const T b,
									 const T c
									 ) noexcept;
	 
	/// \brief Calculate volume of an elliptic cylinder.
	/// Dimensions of the cylinder semi-major axes are \p a and \p b, the height is \p h.
	static constexpr T elliptic_cylinder_vol(const T a, const T b, const T h) noexcept;
	
	/// \brief Calculate volume of an ellipsoidal cap.
	/// Semi-major axes of the ellipsoid are \p a , \p b , and \p c
	/// The cap height is \p h : |h| < c.
	static constexpr T ellipsoid_cap_vol(const T a,
										 const T b,
										 const T c,
										 const T h
										 ) noexcept;
	 
	/// \brief Calculate base area of an ellipsoidal cap.
	/// Dimensions of semi-major axes of the ellipsoid are \p a , \p b , and \p c
	/// The cap height is \p h : |h| < c.
	static constexpr T ellipsoid_cap_base_area(const T a,
											   const T b,
											   const T c,
											   const T h
											   ) noexcept;
	
	/// Calculate surface area of a spheroid.
	/// Spheroid is given by \p r  = {a, b, c}, a = b, i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1.
	static constexpr T spheroid_surf_area(const A3<T>& r,	///< Spheroid dimensions at semi-major axes.
										  Msgr& msgr		///< Printing utility.
										  ) noexcept;
	
	/// \brief Unit normal on surface of an axis-aligned ellipsoid.
	/// Calculate unit normal vector at point \p p on surface
	/// of an axis-aligned ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1 with dimensions \p r = {a,b,c}.
	static constexpr A3<T> unormal_on_ellipsoid(const A3<T>& r,	///< Dimensions of an ellipsoid.
												const A3<T>& p	///< Point on ellipsoid surface.
												) noexcept;

	/// \brief Determine symmetry axes of a spheroid.
	/// The spheroid should be axis-aligned. Returns (-1, -1, -1) if spheroid is a shpere,
	/// otherwise (i, j, k) where (i,j) are axes indexes of unequal dimensions and k is index of the pole axis.
	/// \param r Dimensions of the spheroid semi-major axes.
	static A3<int> spheroid_axes_symmetry(const A3<T>& r) noexcept;

	/// \brief Ellipse resulting from the horizontal plane cross-section of an ellipsoid.
	/// Calculates dimensions {a, b} at semi-axes of an ellipse x^2/a^2 + y^2/b^2 = 1 resulting from
	/// the horizontal z = h plane cross-section of an ellipsoid x^2/e[0]^2 + y^2/e[1]^2 + z^2/e[2]^2 = 1
	/// having dimensions \p e.
	/// \return Semi-axes of the ellipsoid cross-section.
	static constexpr A2<T> ellipsoid_horizontal_crosection0(const A3<T>& e,	///< Dimensions of the ellipsoid semi-axes.
															const T& h		///< z-coordinate of the horizontal plane.
															) noexcept;


	// Line intersections ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// \brief Intersection of a line and an ellipsoid.
	/// Find intersection of a line through a point \p v in the direction unit vector \p d and an ellipsoid.
	// The ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 is given by its semi-axes e = {a, b, c}.
	static constexpr T intersection_line_ellipsoid(const A3<T>& v,	///< Point on a line.
												   const A3<T>& e,	///< Ellipsoid semi-axes.
												   const A3<T>& d	///< Direction of the line.
												   ) noexcept;

	/// \brief Intersection of a line and an ellipse.
	/// Find intersection of a line and a not rotated ellipse centered at the origin.
	/// The line is given by a point \p v in plane and a direction vector \p d = (d0, d1)
	/// The ellipse has dimensions \p e = (a, b): x^2 / a^2 + y^2 / b^2 = 1.
	static constexpr T intersection_line_ellipse(const A2<T>& v,	///< Point on the line.
												 const A2<T>& e,	///< Dimensions of the ellipse semi-major axes.
												 const A2<T>& d		///< Line direction vector.
												 ) noexcept;

	/// \brief Intersection of a line and a rotated ellipse centered at the origin.
	/// The line is given by a point \p v = (v0, v1) in plane and a direction vector \p d = (d0, d1).
	/// The ellipse is rotated counterclockwise through angle alpha about the origin, has semi-axes \p e = (a, b)
	/// (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1.
	/// \see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
	static constexpr T intersection_line_ellipse(const A2<T>& v,	///< Point on the line.
												 const A2<T>& d,	///< Line direction vector.
												 const A2<T>& e,	///< Dimensions of the ellipse semi-major axes.
												 const T alpha
												 ) noexcept;

	// Distance of a point p from a plane given by eq. dot(n,x) + o = 0.
//	static constexpr T distance_point_plane( const A3<T>& p, const A3<T> n, const A3<T> o ) noexcept;
	
	// Intersection of a line segment and a plane.
//	static constexpr bool intersection_segment_plane( const A3<T>& p1, const A3<T>& p2, const A3<T>& n,  const A3<T>& o, A3<T>& intersP ) noexcept;
	
	/// \brief Find intersection of a line and a plane.
	static constexpr T intersection_line_plane(const A3<T>& pab,
											   const A3<T>& p10,
											   const A3<T>& p20,
											   const A3<T>& pa0,
											   const T s,
											   Msgr &msgr
											   ) noexcept;

	/// \brief Find intersection of a line and a plane.
	/// The line is defined by a point \p p0 and direction vector \p d .
	/// The plane is defined by three points \p v1  \p v2  \p v3
	/// \return Distance in direction \p d from \p p0 to the intersection point
	static constexpr T intersection_line_plane(const A3<T>& p0,		///< Point on the line.
											   const A3<T>& d,		///< Line direction vector.
											   const A3<T>& v1,		///< Point on the plane.
											   const A3<T>& v2,		///< Point on the plane.
											   const A3<T>& v3		///< Point on the plane.
											   ) noexcept;

	/// \brief Find intersection of a line and a plane.
	/// The line is defined by a point \p p0 and direction vector \p d
	/// The plane is defined by a point \p v and a normal \p n
	/// \return Distance in direction \p d from \p p0 to the intersection point
	static constexpr T intersection_line_plane(const A3<T>& p0,		///< Point on the line.
											   const A3<T>& d,		///< Line direction vector.
											   const A3<T>& v,		///< Point on the plane.
											   const A3<T>& n		///< Plane unit normal vector.
											   ) noexcept;

	/// \brief Find intersection of a vertical line and a plane.
	/// The line is defined by a point \p p0.
	/// The plane is defined by a point \p v, and a unit normal vector \p n
	/// \return Distance between \p p0 and the intersection point ( parallel or antiparallel to the line depending on \p sign )
	static constexpr T intersection_vertical_line_plane( const A3<T>& p0,	///< Point on the line.
														 const int sign,	///< Directionality (-1, 1) of the result relative to \p d
														 const A3<T>& v,	///< Point on the plane.
														 const A3<T>& n		///< Plane unit normal vector.
														 ) noexcept;

	/// \brief Find intersection of a line and a cone.
	static constexpr T intersection_line_cone( const A3<T>& w,
											   const A3<T>& q,
											   const A3<T>& h,
											   const A3<T>& m,
											   const A3<T>& p,
											   const A3<T>& d
											   ) noexcept;


	// Rotations +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// Calculate rotation matrix for rotation over \p angle around a general axis \p n .
	/// \param[in] n Direction of rotation axis.
	/// \param[in] angle Rotation angle.
	/// \param[out] rm Rotation matrix.
	static constexpr void rotmat(const A3<T> n, const T angle, T rm[3][3]) noexcept;

	/// Calculate rotation matrix for rotation over \p angle around an axis parallel to 'x'.
	/// \param[in] angle Rotation angle.
	/// \param[out] rm Rotation matrix.
	static constexpr void rotmatx(const T angle, T rm[3][3]) noexcept;

	/// Calculate rotation matrix for rotation over \p angle around an axis parallel to 'y'.
	/// \param[in] angle Rotation angle.
	/// \param[out] rm Rotation matrix.
	static constexpr void rotmaty(const T angle, T rm[3][3]) noexcept;

	/// Calculate rotation matrix for rotation over \p angle around an axis parallel to 'z'.
	/// \param[in] angle Rotation angle.
	/// \param[out] rm Rotation matrix.
	static constexpr void rotmatz(const T angle, T rm[3][3]) noexcept;


	// Comparisons +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// Given the dimensions \p e of a 3D body, determine if all sides are equal.
	static constexpr bool all_sides_are_equal(const A3<T>& e) noexcept;
	
	/// Given the dimensions \p e of a 3D body, determine if two sides are equal.
	static constexpr bool two_sides_are_equal(const A3<T>& e) noexcept;
	

	// Some conversions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	///  Convert sphericl coordinates to cartesian coordinates.
	static constexpr A3<T> sph2cart(const T ph,			///< Inclination.
						  			const T th			///< Azimuth.
						  			const T rad=one<T>	///< Radius.
						  			);

	/// Convert spherical coordinates to cartesian coordinates.
	static constexpr A3<T> sphere2cart(const T r,		///< Radius.
									   const T theta,
									   const T sinPhi,
									   const T cosPhi
									   ) noexcept;

	/// Convert polar coordinates to cartesian coordinates.
	static constexpr A2<T> polar2cart(const T r,		///< Radius.
									  const T theta		///< Angle.
									  ) noexcept;

	/// Convert Grad to Rad.
	static constexpr T grad2rad(const T grad) noexcept;

	/// Convert Rad to Grad.
	static constexpr T rad2grad(const T rad) noexcept;


	// Closest points ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// \brief Point on a line closest to the origin.
	/// Find point closest to the origin on the line that passes through point \p p in the direction \p d .
	/// The 3D line is defined with 6 Plücker coordinates L = (d, p × d),
	/// where \p d is the direction of the line, and \p p is any point along the line.
	/// \see https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
	/// \return Point on the line closest to the origin.
	static constexpr A3<T> ptClosest2orgn(const A3<T>& p,	///< Any point along the line.
										  const A3<T>& d	///< Direction of the line.
										  ) noexcept;
	
	/// \brief Point on an ellipsoid closest to line.
	/// Find point closest to line on an ellipsoid given by its semi-major axes \p e
	/// Center of the ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1 with dimensions e = {a,b,c} is at the origin.
	/// https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
	/// \return Point on an ellipsoid closest to line.
	static constexpr A3<T> ellipsoid_closest_point2Line(const A3<T>& ptc2o,	///< point on a 3D line closest to the origin.
														const A3<T>& e		///< Ellipsoid given by its semi-major axes.
														) noexcept;
	

	// Projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// \brief Projection of a vector \p d to z=0 plane.
	static constexpr A3<T> proj2z0plane( const A3<T>& d ) noexcept;

	/// \brief Projection of vector \p v on a plane defined by the normal \p n.
	static constexpr A3<T> vector_proj2plane(const A3<T>& v,
											 const A3<T>& n
											 ) noexcept;
	

	// Hexagonal lattice +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// Find coordinates of a hexagonal lattice centered at \p orig with \p step and number of layers \p numLayers.
	static constexpr std::vector<A2<T>> hexagonal_lattice(const A2<T> orig,
														  const T step,
														  const szt numLayers
														  ) noexcept;

	/// Find number of layers in a hexagonal lattice having \p numVertices vertexes.
	static constexpr szt numLayers_hexagonal_lattice (const szt numVertices) noexcept;


	// Two segments ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// Find cosine of an angle between two segments given by their end points.
	static constexpr T cos_two_segments(const A3<T>& p1,
							 			const A3<T>& p2,
							 			const A3<T>& p3,
							 			const A3<T>& p4
							 			) noexcept;

	/// Find dot product of vectors defined by two segments given by their end points.
	static constexpr T dotpr_two_segments(const A3<T>& p1,
										  const A3<T>& p2,
										  const A3<T>& p3,
										  const A3<T>& p4
										  ) noexcept;

	/// Squared distance between two line segments [\p S10, \p S11] and [\p S20, \p S21] in 3D.
	static constexpr T squared_dist3D_Segment_to_Segment(const A3<T>& S10, const A3<T>& S11,
														 const A3<T>& S20, const A3<T>& S21
														 ) noexcept;

	// Points inside triangle ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	/// \brief Determines if a point is inside a 2D triangle.
	/// Find out if point \p pt is inside triangle given by its vertexes \p v1, \p v2 and \p v3.
	/// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
	static bool point_in_triangle(const A2<T>& pt,
								  const A2<T>& v1,
								  const A2<T>& v2,
								  const A2<T>& v3
								  ) noexcept;

	/// \brief Determines if a point is inside a 2D triangle.
	/// Find out if point \p p is inside triangle given by its vertexes \p v1, \p v2 and \p v3.
	/// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
	static bool point_in_triangle (const T* p,
								   const T* v1,
								   const T* v2,
								   const T* v3
								   ) noexcept;
};

// IMPLEMENTATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

template <typename T> inline
A2<T> Geometric<T>::
ellipse( const T alpha,
		 const A2<T>& ab ) noexcept
{
	return { ab[0]*std::cos(alpha),
			 ab[1]*std::sin(alpha) };
}
//----------------------------------------------------------------------------------------------------------------------

// ellipse area
template <typename T> constexpr
T Geometric<T>::
ellipse_area( const T a,
			  const T b ) noexcept
{
	return pi<T>*a*b;
}

//----------------------------------------------------------------------------------------------------------------------

// ellipsoidal volume
template <typename T> constexpr
T Geometric<T>::
ellipsoid_vol( const T a,
			   const T b,
			   const T c ) noexcept
{
	return four<T>/three<T> * pi<T>*a*b*c;
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
bool Geometric<T>::
is_inside_ellipsoid( const A3<T>& p,
					 const A3<T>& e ) noexcept
{
	return p[0]*p[0]/e[0]/e[0] + p[1]*p[1]/e[1]/e[1] + p[2]*p[2]/e[2]/e[2] < one<T> - EPS<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// elliptic cylinder volume
template <typename T> constexpr
T Geometric<T>::
elliptic_cylinder_vol( const T a,
					   const T b,
					   const T h ) noexcept
{
	return h * ellipseArea(a, b);
}

//----------------------------------------------------------------------------------------------------------------------

// ellipsoidal cap volume for cap height |h| < c
template <typename T> constexpr
T Geometric<T>::
ellipsoid_cap_vol( const T a,
				   const T b,
				   const T c,
				   const T h ) noexcept
{	
	return pi<T> * a*b/(c*c) * h*h * (c - h/three<T>);
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
T Geometric<T>::
ellipsoid_cap_base_area( const T a,
						 const T b,
						 const T c,
						 const T h ) noexcept
{
	// ellipsoidal cap base area for cap height |h| < c 
	return pi<T> * a*b/(c*c) * h * (two<T>*c - h);
}

//----------------------------------------------------------------------------------------------------------------------

// spheroidal surface area
// spheroid is given by r[0:2] = {a, b, c}, a = b, i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1
template <typename T> constexpr
T Geometric<T>::
spheroid_surf_area( const A3<T>& r, Msgr &msgr ) noexcept
{
	if (r[0] == r[1]) {
		const auto a2 = r[0]*r[0]; 
		const auto c2 = r[2]*r[2];
		
		if (r[0] < r[2]) {					// prolate spheroid
			const auto e = std::sqrt(one<T> - a2/c2);
			return twopi<T> * (a2 + std::asin(e) * r[0]*r[2] / e);
		}
		if (r[0] > r[2]) {					// oblate spheroid
			const auto e = std::sqrt( one<T> - c2 / a2 );
			return pi<T> * (two<T>*a2 + std::log((one<T> + e) / (one<T> - e) ) * c2 / e);
		}
		return four<T> * pi<T> * a2;		// sphere
	}

	msgr.print("Error in spheroid_surf_area: spheroid r[0] == r[1] is required", 1);
	std::exit(EXIT_FAILURE);
}

//----------------------------------------------------------------------------------------------------------------------

// determine spheroid axes symmetry from its dimensions r:
// returns (-1, -1, -1) if spheroid is a shpere
// otherwise returns (i, j, k) where (i,j) are axes indexes of unequal dimensions and k is index of the pole axis
template <typename T> inline
A3<int> Geometric<T>::
spheroid_axes_symmetry( const A3<T>& r ) noexcept
{
	if (all_sides_are_equal(r))	return {-1};	// is a sphere

	if (r[0] == r[1]) return {0, 2, 2};			// pole is along 2
	if (r[0] == r[2]) return {0, 1, 1};			// pole is along 1
	if (r[1] == r[2]) return {0, 1, 0};			// pole is along 0
	XASSERT(false, " Geometric::spheroid_axes_symmetry failed");
	return 0;
}

//----------------------------------------------------------------------------------------------------------------------

// returns dimensions {a, b} of an ellipse x^2/a^2 + y^2/b^2 = 1
// resulting from the crossection of an ellipsoid x^2/e[0]^2 + y^2/e[1]^2 + z^2/e[2]^2 = 1 with plane z = h
template <typename T> constexpr
A2<T> Geometric<T>::
ellipsoid_horizontal_crosection0( const A3<T>& e, const T& h ) noexcept
{
	const auto u = std::sqrt(e[2]*e[2] - h*h) / e[2];
	return { u * e[0],
			 u * e[1] };
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and an ellipsoid
// d: the line is given by a point v and a direction unit vector d
// e are the ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 semiaxes e = {a, b, c},
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipsoid( const A3<T>& v,
							 const A3<T>& e,
							 const A3<T>& d ) noexcept
{
	const auto u1 = d[0]*d[0] * e[1]*e[1] * e[2]*e[2] +
					d[1]*d[1] * e[0]*e[0] * e[2]*e[2] +
					d[2]*d[2] * e[0]*e[0] * e[1]*e[1];
	const auto u2 = two<T> * ( d[0] * v[0] * e[1]*e[1] * e[2]*e[2] +
							   d[1] * v[1] * e[0]*e[0] * e[2]*e[2] +
							   d[2] * v[2] * e[0]*e[0] * e[1]*e[1] );
	const auto u3 = v[0]*v[0] * e[1]*e[1] * e[2]*e[2] +
				    v[1]*v[1] * e[0]*e[0] * e[2]*e[2] +
					v[2]*v[2] * e[0]*e[0] * e[1]*e[1] -
					e[0]*e[0] * e[1]*e[1] * e[2]*e[2];
	const auto discr = u2 * u2 - four<T> * u1 * u3;
	
	if (discr >= zero<T>) {										// there is an intersection possible
		const auto t1 = (-u2 - std::sqrt(discr)) / (two<T> * u1);	// putative intersection points
		const auto t2 = (-u2 + std::sqrt(discr)) / (two<T> * u1);
		
		if (t1 > zero<T> && t2 > zero<T>)						// both t1 and t2 are in the positive half-line: intersection is at the closest of t1, t2
			return (t1 <= t2) ? t1 : t2;
		if (t1 > zero<T>) return t1; 							// only t1 is in the positive half-line: intersection is at t1
		if (t2 > zero<T>) return t2;							// only t2 is in the positive half-line: intersection is at t2
		return huge<T>;
	}
	else return huge<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and a not rotated ellipse centered at the origin
// the line is given by a point 'v' = (v0, v1) in plane and a direction vector 'd' = (d0, d1)
// the ellipse has dimensions 'ab' = (a, b): x^2 / a^2 + y^2 / b^2 = 1
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipse( const A2<T>& v,
						   const A2<T>& ab,
						   const A2<T>& d ) noexcept
{
	const auto a2 = ab[0]*ab[0];
	const auto b2 = ab[1]*ab[1];

	const auto u1 = d[0]*d[0] * b2 +
					d[1]*d[1] * a2;
	const auto u2 = two<T> * (d[0] * v[0] * b2 +
							  d[1] * v[1] * a2);
	const auto u3 = v[0]*v[0] * b2 +
					v[1]*v[1] * a2 - a2 * b2;
	const auto discr = u2*u2 - four<T> * u1 * u3;
	
	if (discr >= zero<T>) {								// there is an intersection possible
		const auto sd = std::sqrt(discr);
		const auto t1 = (- u2 - sd) / (two<T> * u1);	// putative intersection points
		const auto t2 = (- u2 + sd) / (two<T> * u1);

		if (t1 > zero<T> &&
			t2 > zero<T>)	// both t1 and t2 are in the positive half-line: intersection is at the closest of t1, t2
			return (t1 <= t2) ? t1 : t2;
		if (t1 > zero<T>) return t1;		// only t1 is in the positive half-line: intersection is at t1
		if (t2 > zero<T>) return t2;		// only t2 is in the positive half-line: intersection is at t2
	}
	return huge<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and a rotated ellipse centered at the origin
// the line is given by a point 'v' = (v0, v1) in plane and a direction vector 'd' = (d0, d1)
// the ellipse is rotated counterclockwise through angle alpha about the origin, has semi-axes 'ab' = (a, b)
// (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1
// see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipse( const A2<T>& v,
						   const A2<T>& d,
						   const A2<T>& ab,
						   const T alpha ) noexcept
{
	const auto sia = std::sin(alpha);
	const auto coa = std::cos(alpha);
	const auto sia2 = sia*sia;
	const auto coa2 = coa*coa;

	const auto a2 = ab[0]*ab[0];
	const auto b2 = ab[1]*ab[1];

	// coefficients of the quadratic form of the ellipse eq: A*x^2 + B*x*y + C*y^2 = 1
	const auto A = coa2/a2 + sia2/b2;
	const auto B = two<T> * coa * sia * (one<T>/a2 - one<T>/b2);
	const auto C = sia2/a2 + coa2/b2;

	// coefficients of the quadratic equation of the line-ellipse intersection
	const auto u1 = d[0]*d[0] * A +
					d[0]*d[1] * B +
					d[1]*d[1] * C;
	const auto u2 = two<T> * (A * v[0] * d[0] +
							  C * v[1] * d[1]) +
					B * (v[0] * d[1] +
						 v[1] * d[0]);
	const auto u3 = A * v[0]*v[0] +
					B * v[0]*v[1] +
					C * v[1]*v[1] - one<T>;
	const auto discr = u2*u2 - four<T> * u1 * u3;

	if (discr >= zero<T>) {								// there is an intersection possible
		const auto sd = std::sqrt(discr);
		const auto t1 = (- u2 - sd) / (two<T> * u1);	// putative intersection points
		const auto t2 = (- u2 + sd) / (two<T> * u1);

		if (t1 > zero<T> &&
			t2 > zero<T>)		// both t1 and t2 are in the positive half-line: intersection is at the closest of t1, t2
			return (t1 <= t2) ? t1 : t2;
		if (t1 > zero<T>) return t1;				// only t1 is in the positive half-line: intersection is at t1
		if (t2 > zero<T>) return t2;				// only t2 is in the positive half-line: intersection is at t2
	}
	return huge<T>;
}

//----------------------------------------------------------------------------------------------------------------------

/*
// Distance of a point p from a plane given by eq. dot(n,x) + o = 0
template <typename T> constexpr
T Geometric<T>::
distance_point_plane( const A3<T>& p, const A3<T> n,  const A3<T> o ) noexcept
{
	return p.dot(n) + o;
}

// intersection of a line segment and and a plane
template <typename T> constexpr
bool Geometric<T>::
intersection_segment_plane( const A3<T>& p1, const A3<T>& p2, const A3<T>& n,  const A3<T>& o, A3<T>& intersP ) noexcept
{
	const auto d1 = distFromPlane(p1, n, o);
	const auto d2 = distFromPlane(p2, n, o);

	const bool bP1OnPlane = (std::abs(d1) < EPS<T>);
	const bool bP2OnPlane = (std::abs(d2) < EPS<T>);

	if (bP1OnPlane) {
		intersP = p1;
		return true;
	}

	if (bP2OnPlane) {
		intersP = p2;
		return true;
	}

	if (bP1OnPlane && bP2OnPlane)
		return false;

	if (std::signbit(d1) == std::signbit(d2))  // points on the same side of plane
		return false;

	const auto t = d1 / (d1 - d2); // 'time' of intersection point on the segment
	intersP = p1 + t * (p2 - p1);
	return true;
}
*/

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and and a plane
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane( const A3<T>& pab,
						 const A3<T>& p10,
						 const A3<T>& p20,
						 const A3<T>& pa0,
						 const T s,
						 Msgr& msgr ) noexcept
{
	const auto den1 = p10[1]*p20[0] -
					  p20[1]*p10[0];
	
	const auto kappa = (p20[1]*pab[0] - pab[1]*p20[0]) / den1;
	const auto lamda = (pa0[1]*p20[0] - p20[1]*pa0[0]) / den1;
	
	const auto den2 = p20[0] * (pab[2] + p10[2]*kappa) -
					  p20[2] * (pab[0] + p10[0]*kappa);
	
	if (!den1 || !den2)
		msgr.exit("Error: singularity in line-plane: filament");

	return ( p20[0] * (pa0[2] - p10[2]*lamda) -
			 p20[2] * (pa0[0] - p10[0]*lamda) ) / den2 * s;
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and and a plane
// the line is defined by a point 'p' and direction vector 'd'
// the plane is defined by three points 'v1', 'v2', 'v3'
// returns distance in direction 'd' from 'p0' to the intersection point
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane( const A3<T>& p,
						 const A3<T>& d,
						 const A3<T>& v1,
						 const A3<T>& v2,
						 const A3<T>& v3) noexcept
{
	const auto v13 = v1 - v3;
	const auto v23 = v2 - v3;

	const auto vn = A3<T>::crosspr(v13, v23);
	const auto n = vn.unitv();
	const auto nd = A3<T>::dotpr(n, d);

	if (std::abs(nd) < EPS<T>)   // the line is parallel to the plane
		return -huge<T>;

	const auto nw = -A3<T>::dotpr(n, p - v1);

	return nw / nd;	// intersection point is at: p + t*d
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and and a plane
// the line is defined by a point 'p' and direction vector 'd'
// the plane is defined by a points 'v', and a normal 'n'
// returns distance in direction 'd' from 'p0' to the intersection point
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane( const A3<T>& p,
						 const A3<T>& d,
						 const A3<T>& v,
						 const A3<T>& n) noexcept
{
	const auto nd = A3<T>::dotpr(n, d);

	if (std::abs(nd) < EPS<T>)   // the line is parallel to the plane
		return -huge<T>;

	return -A3<T>::dotpr(n, p - v) / nd;	// intersection point is at: p + t*d
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a vertical line and and a plane
// the line is defined by a point 'p' and direction vector 'd'
// the plane is defined by a points 'v', and a normal 'n'
// returns distance in direction 'd' from 'p0' to the intersection point
template <typename T> constexpr
T Geometric<T>::
intersection_vertical_line_plane( const A3<T>& p,
						 		  const int sign,
						 		  const A3<T>& v,
						 		  const A3<T>& n) noexcept
{
	if (std::abs(n[2]) < EPS<T>)   // the line is parallel to the plane
		return -huge<T>;

	return -A3<T>::dotpr(n, p - v) / sign*n[2];	// intersection point is at: p + t*d
}

//----------------------------------------------------------------------------------------------------------------------

// intersection of a line and and a cone
template <typename T> constexpr
T Geometric<T>::
intersection_line_cone( const A3<T>& w,
						const A3<T>& q,
						const A3<T>& h,
						const A3<T>& m,
						const A3<T>& p,
						const A3<T>& d ) noexcept
{
	const auto u1 = h.dotpr(d*d) + two<T>*(q[0]*q[1]*d[0]*d[1] +
										   q[0]*q[2]*d[0]*d[2] +
									 	   q[1]*q[2]*d[1]*d[2]);
									 
	const auto u2 = m.dotpr(d) + two<T>*(h.dotpr(d*p) + q[0]*q[1]*(d[0]*p[1]+d[1]*p[0]) +
													    q[0]*q[2]*(d[0]*p[2]+d[2]*p[0]) +
													    q[1]*q[2]*(d[1]*p[2]+d[2]*p[1]));
													
	const auto u3 = m.dotpr(p) + h.dotpr(p*p+w*w) + two<T>*(q[0]*q[1]*(p[0]*p[1]+w[0]*w[1]) +
															q[0]*q[2]*(p[0]*p[2]+w[0]*w[2]) +
															q[1]*q[2]*(p[1]*p[2]+w[1]*w[2]));
	auto discr = u2 * u2 - four<T> * u1 * u3;
	
	if (discr >= zero<T>) {											// there is an intersection possible
		const auto t1 = (- u2 - std::sqrt(discr)) / (two<T> * u1);	// putative intersection points
		const auto t2 = (- u2 + std::sqrt(discr)) / (two<T> * u1);

		if(		 t1 > zero<T> && 
				 t2 > zero<T>) return (t1 <= t2) ? t1 : t2;	// both t1 and t2 are in the positive half-line: intersection is at the closest of t1, t2
		else if (t1 > zero<T>) return t1; 					// only t1 is in the positive half-line: intersection is at t1
		else if (t2 > zero<T>) return t2;					// only t2 is in the positive half-line: intersection is at t2
		else				   return -one<T>;
	}
	else return -one<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// projection of a vector 'd' to z=0 plane
template <typename T> constexpr inline
A3<T> Geometric<T>::
proj2z0plane( const A3<T>& d ) noexcept
{
	return {d(0,1).unitv(), zero<T>};
}

//----------------------------------------------------------------------------------------------------------------------

// returns as 'rm' a rotation matrix for rotation over 'angle' around a general axis 'n'
template <typename T> constexpr
void Geometric<T>::
rotmat( const A3<T> n, T angle, T rm[3][3] ) noexcept
{
	const auto sia = std::sin(angle);
	const auto coa = std::cos(angle);
	const auto c = - n*n + one<T>;
	
	rm[0][0] = n[0]*n[0] + c[0]*coa;				rm[0][1] = n[0]*n[1]*(one<T>-coa) - n[2]*sia;	rm[0][2] = n[0]*n[2]*(one<T>-coa) + n[1]*sia;
	rm[1][0] = n[0]*n[1]*(one<T>-coa) + n[2]*sia;	rm[1][1] = n[1]*n[1] + c[1]*coa;				rm[1][2] = n[1]*n[2]*(one<T>-coa) - n[0]*sia;
	rm[2][0] = n[0]*n[2]*(one<T>-coa) - n[1]*sia;	rm[2][1] = n[1]*n[2]*(one<T>-coa) + n[0]*sia;	rm[2][2] = n[2]*n[2] + c[2]*coa;
}

//----------------------------------------------------------------------------------------------------------------------

// returns as 'rm' a rotation matrix for rotation over 'angle' around an axis parallel to 'x'
template <typename T> constexpr
void Geometric<T>::
rotmatx( T angle, T rm[3][3] ) noexcept
{
	const A3<T> n {one<T>, zero<T>, zero<T>};
	
	const auto sia = std::sin(angle);
	const auto coa = std::cos(angle);
	
	rm[0][0] = n[0]*n[0] + (one<T>-n[0]*n[0])*coa;	rm[0][1] = zero<T>;		rm[0][2] = zero<T>; 
	rm[1][0] = zero<T>;								rm[1][1] = zero<T>;		rm[1][2] = - n[0]*sia;
	rm[2][0] = zero<T>;								rm[2][1] = n[0]*sia;	rm[2][2] = zero<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// returns as 'rm' a rotation matrix for rotation over 'angle' around an axis parallel to 'y'
template <typename T> constexpr
void Geometric<T>::
rotmaty( T angle, T rm[3][3] ) noexcept
{
	const A3<T> n {zero<T>, one<T>, zero<T>};
	
	const auto sia = std::sin(angle);
	const auto coa = std::cos(angle)}
	
	rm[0][0] = zero<T>;		rm[0][1] = zero<T>;								rm[0][2] = n[1]*sia; 
	rm[1][0] = zero<T>;		rm[1][1] = n[1]*n[1] + (one<T>-n[1]*n[1])*coa;	rm[1][2] = zero<T>;
	rm[2][0] = - n[1]*sia;	rm[2][1] = zero<T>;								rm[2][2] = zero<T>;
}

//----------------------------------------------------------------------------------------------------------------------

// returns as 'rm' a rotation matrix for rotation over 'angle' around an axis parallel to 'z'
template <typename T> constexpr
void Geometric<T>::
rotmatz( T angle, T rm[3][3] ) noexcept
{
	const A3<T> n {zero<T>, zero<T>, one<T>};
	
	const auto sia = std::sin(angle);
	const auto coa = std::cos(angle);
	
	rm[0][0] = zero<T>;		rm[0][1] = - n[2]*sia;	rm[0][2] = zero<T>; 
	rm[1][0] = n[2]*sia;	rm[1][1] = zero<T>;		rm[1][2] = zero<T>;
	rm[2][0] = zero<T>;		rm[2][1] = zero<T>;		rm[2][2] = n[2]*n[2] + (one<T>-n[2]*n[2])*coa;
}

//----------------------------------------------------------------------------------------------------------------------

// unit normal at point p on surface of an axis-aligned ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1 with dimensions r = {a,b,c}
template <typename T> constexpr
A3<T> Geometric<T>::
unormal_on_ellipsoid( const A3<T>& r, const A3<T>& p ) noexcept
{	
	const auto a2 = r[0]*r[0];
	const auto b2 = r[1]*r[1];
	const auto c2 = r[2]*r[2];
	const A3<T> n {b2*c2*p[0],
				   a2*c2*p[1],
				   a2*b2*p[2]};
	return n.unitv();
}

//----------------------------------------------------------------------------------------------------------------------

// given the dimensions of a 3D body, determine if all_sides_are_equal
template <typename T> constexpr inline
bool Geometric<T>::
all_sides_are_equal( const A3<T>& e ) noexcept
{
	return e[0] == e[1] &&
		   e[1] == e[2];
}

//----------------------------------------------------------------------------------------------------------------------

// given the dimensions of a 3D body, determine if two_sides_are_equal
template <typename T> constexpr inline
bool Geometric<T>::
two_sides_are_equal( const A3<T>& e ) noexcept
{
	return e[0] == e[1] ||
		   e[0] == e[2] ||
		   e[1] == e[2];
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> inline
A3<T> Geometric<T>::
sph2cart( const T ph,		// inclination
		  const T th,		// azimuth
		  const T rad )
{
	return { std::cos(ph) * std::cos(th),
			 std::cos(ph) * std::sin(th),
			 std::sin(ph) } * rad;
}

//----------------------------------------------------------------------------------------------------------------------

// conversion of spherical to cartesian coordinates
// phi: inclination;   theta: azimuth
template <typename T> constexpr inline
A3<T> Geometric<T>::
sphere2cart( const T r,
			 const T theta,
			 const T sinPhi,
			 const T cosPhi ) noexcept
{
	return { r * std::sin(theta) * sinPhi,
			 r * std::cos(theta) * sinPhi,
			 r * cosPhi };
}

//----------------------------------------------------------------------------------------------------------------------

// conversion of polar to cartesian coordinates
// phi: inclination;   theta: azimuth
template <typename T> constexpr inline
A2<T> Geometric<T>::
polar2cart( const T r,
			const T theta ) noexcept
{
	return { r * std::sin(theta),
			 r * std::cos(theta) };
}

//----------------------------------------------------------------------------------------------------------------------

// point on the line closest to the origin
template <typename T> constexpr inline
A3<T> Geometric<T>::
ptClosest2orgn( const A3<T>& p,
				const A3<T>& d ) noexcept
{
	// https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
	// A 3D line is defined with 6 Plücker coordinates L = ( d, p × d ) where d is the direction of the line, and p is any point along the line

	return A3<T>::crosspr(d, A3<T>::crosspr(d, p));
}

//----------------------------------------------------------------------------------------------------------------------

// point on an ellipsoid closest to line
template <typename T> constexpr inline
A3<T> Geometric<T>::
ellipsoid_closest_point2Line( const A3<T>& ptc2o,
							  const A3<T>& elps ) noexcept
{
	// https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
	// ptc2o is a point on a 3D line closest to the origin; 
	// center of the ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1 with dimensions elps = {a,b,c} is at the origin

	const auto elps2 = elps * elps;
	const auto l = std::sqrt(A3<T>::dotpr(ptc2o*ptc2o, elps2));

	return (ptc2o * elps2) / l;
}

//----------------------------------------------------------------------------------------------------------------------

// projection of vector v on a plane given by the normal n
template <typename T> constexpr inline
A3<T> Geometric<T>::
vector_proj2plane( const A3<T>& v,
				   const A3<T>& n ) noexcept
{
	return v - v.vecProjection(n);
}

//----------------------------------------------------------------------------------------------------------------------

// squared distance between two line segments in 3D
template <typename T> constexpr
T Geometric<T>::
squared_dist3D_Segment_to_Segment( const A3<T>& S10,
								   const A3<T>& S11,
								   const A3<T>& S20,
								   const A3<T>& S21 ) noexcept
{
	constexpr auto SMALL_NUM = static_cast<T>(0.00000001);		// anything that avoids division overflow
	const auto u = S11 - S10;
	const auto v = S21 - S20;
	const auto w = S10 - S20;
	const auto a = A3<T>::dotpr(u,u);			// always >= 0
	const auto b = A3<T>::dotpr(u,v);
	const auto c = A3<T>::dotpr(v,v);			// always >= 0
	const auto d = A3<T>::dotpr(u,w);
	const auto e = A3<T>::dotpr(v,w);
	const auto D = a*c - b*b;					// always >= 0
	T sc, sN, sD = D;			// sc = sN / sD, default sD = D >= 0
	T tc, tN, tD = D;			// tc = tN / tD, default tD = D >= 0
	
	// compute the line Config of the two closest points
	if (D < SMALL_NUM) {			// the lines are almost parallel
		sN = zero<T>;				// force using point P0 on segment S1
		sD = one<T>;				// to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {							// get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < zero<T> ) {		// sc < 0 => the s=0 edge is visible
			sN = zero<T>;
			tN = e;
			tD = c;
		}
		else if (sN > sD ) {		// sc > 1  => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}
	if (tN < zero<T>) {				// tc < 0 => the t=0 edge is visible
		tN = zero<T>;
		// recompute sc for this edge
		if (-d < zero<T>)	sN = zero<T>;
		else if (-d > a)	sN = sD;
		else {				sN = -d;
							sD = a;
		}
	}
	else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < zero<T>) sN = zero<T>;
		else if ((-d + b) > a)	sN = sD;
		else {					sN = (-d +  b);
								sD = a;
		}
	}
	// finally do the division to get sc and tc
	sc = (std::abs(sN) < SMALL_NUM ? zero<T> : sN / sD);
	tc = (std::abs(tN) < SMALL_NUM ? zero<T> : tN / tD);
	
	// get the difference of the two closest points
	const auto dP = w + (u * sc) - (v * tc);  // =  S1(sc) - S2(tc)
	
	return dP.dotpr();   // return the closest distance
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T>
std::vector<A2<T>> Geometric<T>::
hexagonal_lattice( const A2<T> orig,
				   const T step,
				   const szt numLayers ) noexcept
{
	std::vector<A2<T>> v;

	auto add_point = [&](const szt j, const int sign)
	{
    for (szt i=1; i<=2*numLayers+1-j; i++)
		v.emplace_back(orig + A2<T> {(half<T>*j + i - numLayers - one<T>),
									 sign * static_cast<T>(j) * std::sin(pi<T>/three<T>)} * step);
	};

	point(0, 0);
	for (szt j=1; j<=numLayers; j++) {
		add_point(j, 1);
		add_point(j, -1);
	}
	return v;
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
szt Geometric<T>::
numLayers_hexagonal_lattice ( const szt numVertices ) noexcept
{
	if (numVertices <= 1)
		return numVertices;

	szt numLayers {1};
	szt n {1};
	do
		n += 6 * numLayers++;
	while (n < numVertices);

	return numLayers;
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
T Geometric<T>::
grad2rad( const T grad ) noexcept
{
	return grad*pi<T>/static_cast<T>(180.);
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
T Geometric<T>::
rad2grad( const T rad ) noexcept
{
	return rad*static_cast<T>(180.)/pi<T>;
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
T Geometric<T>::
cos_two_segments( const A3<T>& p1,
				  const A3<T>& p2,
				  const A3<T>& p3,
				  const A3<T>& p4 ) noexcept
{
	const auto d1 = p2 - p1;
	const auto d2 = p4 - p3;

	return d1.dotpr(d2) / (d1.norm() * d2.norm());
}

//----------------------------------------------------------------------------------------------------------------------

template <typename T> constexpr
T Geometric<T>::
dotpr_two_segments( const A3<T>& p1,
					const A3<T>& p2,
					const A3<T>& p3,
					const A3<T>& p4 ) noexcept
{
	const auto d1 = p2 - p1;
	const auto d2 = p4 - p3;

	return d1.dotpr(d2);
}

//----------------------------------------------------------------------------------------------------------------------

// Determines if a point is inside a 2D triangle
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template <typename T> inline
bool Geometric<T>::
point_in_triangle (const A2<T>& p, const A2<T>& v1, const A2<T>& v2, const A2<T>& v3) noexcept
{
	auto sign = []( const A2<T>& p1, const A2<T>& p2, const A2<T>& p3 ) noexcept
	{
		const auto d {(p1[0] - p3[0]) * (p2[1] - p3[1]) -
					  (p2[0] - p3[0]) * (p1[1] - p3[1])};
		return d;
	};

    const T d1 {sign(p, v1, v2)};
    const T d2 {sign(p, v2, v3)};
    const T d3 {sign(p, v3, v1)};

    return !(((d1 < zero<T>) || (d2 < zero<T>) || (d3 < zero<T>)) &&
			 ((d1 > zero<T>) || (d2 > zero<T>) || (d3 > zero<T>)));
}

//----------------------------------------------------------------------------------------------------------------------

// Determines if a point is inside a 2D triangle
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template <typename T> inline
bool Geometric<T>::
point_in_triangle (const T* p, const T* v1, const T* v2, const T* v3) noexcept
{
	auto sign = []( const T* p1, const T* p2, const T* p3 ) noexcept
	{
		const auto d {(*p1 - *p3) * (*(p2+1) - *(p3+1)) -
					  (*p2 - *p3) * (*(p1+1) - *(p3+1))};
		return d;
	};

    const T d1 {sign(p, v1, v2)};
    const T d2 {sign(p, v2, v3)};
    const T d3 {sign(p, v3, v1)};

    return !(((d1 < zero<T>) || (d2 < zero<T>) || (d3 < zero<T>)) &&
			 ((d1 > zero<T>) || (d2 > zero<T>) || (d3 > zero<T>)));
}

}	// namespace Common
}	// namespace Utils

#endif 	// UTILS_COMMON_GEOMETRIC_FUNCTIONS
