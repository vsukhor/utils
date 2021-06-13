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
 * \file misc.h
 * \brief A collection of genetally useful geometric functions.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_COMMON_GEOMETRIC_FUNCTIONS
#define UTILS_COMMON_GEOMETRIC_FUNCTIONS

#include "../arrays/all.h"
#include "../msgr.h"
#include "exceptions.h"
#include "misc.h"


/// General stuff.
namespace utils::common {

/**
 * \class Geometric geometric_functions.h
 * \brief A loose collection of geometry-related static functions.
 * \tparam T Floating point type.
 */
template <typename T>
struct Geometric {

    // Make sure that the template parameter is a floating type.
    static_assert(std::is_floating_point<T>::value,
            "Class Geometric can only be instantiated with floating point types");


    using A2t = arrays::A2<T>;
    using A3t = arrays::A3<T>;
    using A3i = arrays::A3<int>;
    using SimpleException = exceptions::Simple;

    static constexpr auto RAD2GRAD = static_cast<T>(180);

    // Elliptic shapes +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Get a point on an ellipse centered at zero.
    static A2t ellipse(
        T alpha,     ///< Angular coordinate.
        const A2t& ab    ///< Dimensions of ellipse semi-major axes.
    ) noexcept;

    /// \brief Point inside an ellipsoid or outside.
    /// \details Find out if point \p p is inside an ellipsoid given by
    /// dimensions \p e of the semi-major axes.
    static constexpr bool is_inside_ellipsoid(
        const A3t& p,
        const A3t& e
    ) noexcept;
    
    /// \brief Calculate area of an ellipse.
    /// \param a  Length of the semi-major axis.
    /// \param b  Length of the semi-major axis.
    static constexpr T ellipse_area(T a, T b) noexcept;

    /// \brief Calculate volume of an ellipsoid.
    /// \param a  Length of the semi-major axis.
    /// \param b  Length of the semi-major axis.
    /// \param c  Length of the semi-major axis.

    static constexpr T ellipsoid_vol(T a, T b, T c) noexcept;
     
    /// \brief Calculate volume of an elliptic cylinder.
    /// \param a  Length of the semi-major axis.
    /// \param b  Length of the semi-major axis.
    /// \param h Cylinder height.
    static constexpr T elliptic_cylinder_vol(T a, T b, T h) noexcept;
    
    /// \brief Calculate volume of an ellipsoidal cap.
    /// \param a  Length of the semi-major axis.
    /// \param b  Length of the semi-major axis.
    /// \param c  Length of the semi-major axis.
    /// \param h The cap height (|h| < c).
    static constexpr T ellipsoid_cap_vol(T a, T b, T c, T h) noexcept;
     
    /// \brief Calculate base area of an ellipsoidal cap.
    /// \param a  Length of the semi-major axis.
    /// \param b  Length of the semi-major axis.
    /// \param c  Length of the semi-major axis.
    /// \param h The cap height (|h| < c).
    static constexpr T ellipsoid_cap_base_area(T a, T b,  T c, T h) noexcept;
    
    /// \brief Calculate surface area of a spheroid.
    /// \details Spheroid is given by \p
    /// r  = {a, b, c}, a = b, i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1.
    static constexpr T spheroid_surf_area(
        const A3t& r,    ///< Spheroid dimensions at semi-major axes.
        Msgr& msgr       ///< Output message processor.
    );

    /// \brief Unit normal on the surface of an axis-aligned ellipsoid.
    /// \details Calculate unit normal vector at point \a p on surface
    /// of an axis-aligned ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    /// with dimensions \p r = {a,b,c}.
    static constexpr auto unormal_on_ellipsoid(
        const A3t& r,    ///< Dimensions of an ellipsoid.
        const A3t& p     ///< Point on ellipsoid surface.
        ) noexcept -> A3t;

    /// \brief Determine symmetry axes of a spheroid.
    /// \details The spheroid should be axis-aligned.
    /// \param r Dimensions of the spheroid semi-major axes.
    /// \return (-1, -1, -1) if spheroid is a shpere,
    /// otherwise (i, j, k) where (i,j) are axes indexes of unequal dimensions
    /// and k is index of the pole axis.
    static auto spheroid_axes_symmetry(
        const A3t& r
    ) noexcept -> A3i;

    /// \brief Ellipse resulting from the horizontal plane cross-section of an ellipsoid.
    /// \details Calculates dimensions {a, b} at semi-axes of an
    /// ellipse x^2/a^2 + y^2/b^2 = 1 resulting from the horizontal z = h plane
    /// cross-section of an ellipsoid x^2/e[0]^2 + y^2/e[1]^2 + z^2/e[2]^2 = 1
    /// having dimensions \p e.
    /// \return Semi-axes of the ellipsoid cross-section.
    static constexpr auto ellipsoid_horizontal_crosection0(
        const A3t& e,    ///< Dimensions of the ellipsoid semi-axes.
        T h              ///< z-coordinate of the horizontal plane.
    ) noexcept -> A2t;


    // Line intersections ++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// \brief Intersection of a line and an ellipsoid.
    /// \details Find intersection of a line through a point \p v in the direction
    /// unit vector \p d and an ellipsoid. The ellipsoid
    /// x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 is given by its semi-axes e = {a, b, c}.
    static constexpr T intersection_line_ellipsoid(
        const A3t& v,    ///< Point on a line.
        const A3t& e,    ///< Ellipsoid semi-axes.
        const A3t& d     ///< Direction of the line.
    ) noexcept;

    /// \brief Intersection of a line and an ellipse.
    /// Find intersection of a line and a not rotated ellipse centered at
    /// the origin. The line is given by a point \p v in plane and a
    /// direction vector \p d = (d0, d1) The ellipse has dimensions
    /// \p e = (a, b): x^2 / a^2 + y^2 / b^2 = 1.
    static constexpr T intersection_line_ellipse(
        const A2t& v,     ///< Point on the line.
        const A2t& ab,    ///< Dimensions of the ellipse semi-major axes.
        const A2t& d      ///< Line direction vector.
    ) noexcept;

    /// \brief Intersection of a line and a rotated ellipse centered at the origin.
    /// \details The line is given by a point \p v = (v0, v1) in plane and a
    /// direction vector \p d = (d0, d1). The ellipse is rotated counterclockwise
    /// through angle alpha about the origin, has semi-axes \p ab = (a, b)
    /// (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1.
    /// \see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
    static constexpr T intersection_line_ellipse(
        const A2t& v,    ///< Point on the line.
        const A2t& d,    ///< Line direction vector.
        const A2t& ab,   ///< Dimensions of the ellipse semi-major axes.
        T alpha
    ) noexcept;

    // Distance of a point p from a plane given by eq. dot(n,x) + o = 0.
    // static constexpr T distance_point_plane(
    //      const A3t& p, const A3t n, const A3t o ) noexcept;

    // Intersection of a line segment and a plane.
    // static constexpr bool intersection_segment_plane(
    //      const A3t& p1, const A3t& p2, const A3t& n,
    //      const A3t& o, A3t& intersP ) noexcept;
    
    /// \brief Find intersection of a line and a plane.
    static constexpr T intersection_line_plane(
        const A3t& pab,
        const A3t& p10,
        const A3t& p20,
        const A3t& pa0,
        T s,
        Msgr &msgr
    ) noexcept;

    /// \brief Find intersection of a line and a plane.
    /// \details The line is defined by a point \p p0 and direction vector \p d .
    /// The plane is defined by three points \p v1  \p v2  \p v3
    /// \return Distance in direction \p d from \p p0 to the intersection point
    static constexpr T intersection_line_plane(
        const A3t& p0,       ///< Point on the line.
        const A3t& d,        ///< Line direction vector.
        const A3t& v1,       ///< Point on the plane.
        const A3t& v2,       ///< Point on the plane.
        const A3t& v3        ///< Point on the plane.
    ) noexcept;

    /// \brief Find intersection of a line and a plane.
    /// \details The line is defined by a point \p p0 and direction vector \p d
    /// The plane is defined by a point \p v and a normal \p n
    /// \return Distance in direction \p d from \p p0 to the intersection point.
    static constexpr T intersection_line_plane(
        const A3t& p0,       ///< Point on the line.
        const A3t& d,        ///< Line direction vector.
        const A3t& v,        ///< Point on the plane.
        const A3t& n         ///< Plane unit normal vector.
    ) noexcept;

    /// \brief Find intersection of a vertical line and a plane.
    /// \details The line is defined by a point \p p0.
    /// The plane is defined by a point \p v, and a unit normal vector \p n
    /// \return Distance between \p p0 and the intersection point (parallel or
    /// antiparallel to the line depending on \p sign)
    static constexpr T intersection_vertical_line_plane(
        const A3t& p0,   ///< Point on the line.
        int sign,        ///< Directionality (-1, 1) of the result relative to \p d
        const A3t& v,    ///< Point on the plane.
        const A3t& n     ///< Plane unit normal vector.
    ) noexcept;

    /// \brief Find intersection of a line and a cone.
    static constexpr T intersection_line_cone(
        const A3t& w,
        const A3t& q,
        const A3t& h,
        const A3t& m,
        const A3t& p,
        const A3t& d
    ) noexcept;


    // Rotations +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// \brief Rotation matrix for rotation over \p angle around a general axis \p n .
    /// \param[in] n Direction of rotation axis.
    /// \param[in] angle Rotation angle.
    /// \param[out] rm Rotation matrix.
    static constexpr void rotmat(
        const A3t& n,
        T angle,
        T rm[3][3]
    ) noexcept;

    /// \brief Rotation matrix for rotation over \p angle around an axis parallel to 'x'.
    /// \param[in] angle Rotation angle.
    /// \param[out] rm Rotation matrix.
    static constexpr void rotmatx(
        T angle,
        T rm[3][3]
    ) noexcept;

    /// \brief Rotation matrix for rotation over \p angle around an axis parallel to 'y'.
    /// \param[in] angle Rotation angle.
    /// \param[out] rm Rotation matrix.
    static constexpr void rotmaty(
        T angle,
        T rm[3][3]
    ) noexcept;

    /// \brief Rotation matrix for rotation over \p angle around an axis parallel to 'z'.
    /// \param[in] angle Rotation angle.
    /// \param[out] rm Rotation matrix.
    static constexpr void rotmatz(
        T angle,
        T rm[3][3]
    ) noexcept;


    // Comparisons +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Given the dimensions \p e of a 3D body, determine if all sides are equal.
    static constexpr bool all_sides_are_equal(const A3t& e) noexcept;
    
    /// Given the dimensions \p e of a 3D body, determine if two sides are equal.
    static constexpr bool two_sides_are_equal(const A3t& e) noexcept;
    

    // Some conversions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Convert spherical coordinates to cartesian coordinates.
    static constexpr auto sph2cart(
        T ph,            ///< Inclination.
        T th,            ///< Azimuth.
        T rad=one<T>     ///< Radius.
    ) -> A3t;

    /// Convert spherical coordinates to cartesian coordinates.
    static constexpr auto sphere2cart(
        T r,        ///< Radius.
        T theta,
        T sinPhi,
        T cosPhi
    ) noexcept  -> A3t;

    /// Convert polar coordinates to cartesian coordinates.
    static constexpr auto polar2cart(
        T r,        ///< Radius.
        T theta     ///< Angle.
    ) noexcept -> A2t;

    /// Convert Grad to Rad.
    static constexpr T grad2rad(T grad) noexcept;

    /// Convert Rad to Grad.
    static constexpr T rad2grad(T rad) noexcept;


    // Closest points ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Point on a line closest to the origin.
    /// Find point closest to the origin on the line that passes
    /// through point \p p in the direction \p d . The 3D line is
    /// defined with 6 Plücker coordinates L = (d, p × d),
    /// where \p d is the direction of the line, and \p p is any point along the line.
    /// \see https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    /// \return Point on the line closest to the origin.
    static constexpr auto ptClosest2orgn(
        const A3t& p,    ///< Any point along the line.
        const A3t& d    ///< Direction of the line.
    ) noexcept -> A3t;
    
    /// \brief Point on an ellipsoid closest to line.
    /// \details Find point closest to line on an ellipsoid given by its
    /// semi-major axes \p e Center of the ellipsoid
    ///      (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    /// with dimensions e = {a,b,c} is at the origin.
    /// https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    /// \return Point on an ellipsoid closest to line.
    static constexpr auto ellipsoid_closest_point2Line(
        const A3t& ptc2o,    ///< point on a 3D line closest to the origin.
        const A3t& e         ///< Ellipsoid given by its semi-major axes.
    ) noexcept  -> A3t;
    

    // Projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Projection of a vector \b d to z=0 plane.
    static constexpr auto proj2z0plane( const A3t& d ) noexcept -> A2t;

    /// Projection of vector \b v on a plane defined by the normal \b n.
    static constexpr auto vector_proj2plane(
        const A3t& v,
        const A3t& n
    ) noexcept -> A3t;
    

    // Hexagonal lattice +++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// \brief Find coordinates of a hexagonal lattice.
    /// \details Find coordinates of a hexagonal lattice centered at \p orig
    /// with \p step and number of layers \p numLayers.
    static constexpr auto hexagonal_lattice(
        A2t orig,
        T step,
        szt numLayers
    ) noexcept -> std::vector<A2t>;

    /// Find number of layers in a hexagonal lattice having \p numVertices vertexes.
    static constexpr auto numLayers_hexagonal_lattice(
        szt numVertices
    ) noexcept -> szt;


    // Two segments ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /// Find cosine of an angle between two segments given by their end points.
    static constexpr T cos_two_segments(
        const A3t& p1,
        const A3t& p2,
        const A3t& p3,
        const A3t& p4
    ) noexcept;

    /// Find dot product of vectors defined by two segments given by their end points.
    static constexpr T dotpr_two_segments(
        const A3t& p1,
        const A3t& p2,
        const A3t& p3,
        const A3t& p4
    ) noexcept;

    /// Squared distance between two line segments [\p S10, \p S11] and [\p S20, \p S21] in 3D.
    static constexpr T squared_dist3D_Segment_to_Segment(
        const A3t& S10, const A3t& S11,
        const A3t& S20, const A3t& S21
    ) noexcept;

    // Points inside triangle ++++++++++++++++++++++++++++++++++++++++++++++++++

    /// \brief Determines if a point is inside a 2D triangle.
    /// \details Find out if point \p pt is inside triangle given by its
    /// vertexes \p v1, \p v2 and \p v3.
    /// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
    static bool point_in_triangle(
        const A2t& pt,
        const A2t& v1,
        const A2t& v2,
        const A2t& v3
    ) noexcept;

    /// \brief Determines if a point is inside a 2D triangle.
    /// \details Find out if point \p p is inside triangle given by its
    /// vertexes \p v1, \p v2 and \p v3.
    /// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
    static bool point_in_triangle (
        const T* p,
        const T* v1,
        const T* v2,
        const T* v3
    ) noexcept;
};


// IMPLEMENTATION xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


template <typename T> inline
auto Geometric<T>::
ellipse(
    const T alpha,
    const A2t& ab
) noexcept -> A2t
{
    return { ab[0]*std::cos(alpha),
             ab[1]*std::sin(alpha) };
}


// ellipse area
template <typename T> constexpr
T Geometric<T>::
ellipse_area(
    const T a,
    const T b
) noexcept
{
    return pi<T>*a*b;
}


// ellipsoidal volume
template <typename T> constexpr
T Geometric<T>::
ellipsoid_vol(
    const T a,
    const T b,
    const T c
) noexcept
{
    return four<T>/three<T> * pi<T>*a*b*c;
}


template <typename T> constexpr
bool Geometric<T>::
is_inside_ellipsoid(
    const A3t& p,
    const A3t& e
) noexcept
{
    return p[0]*p[0]/e[0]/e[0] +
           p[1]*p[1]/e[1]/e[1] +
           p[2]*p[2]/e[2]/e[2] < one<T> - EPS<T>;
}


// Elliptic cylinder volume.
template <typename T> constexpr
T Geometric<T>::
elliptic_cylinder_vol(
    const T a,
    const T b,
    const T h
) noexcept
{
    return h * ellipseArea(a, b);
}


// Ellipsoidal cap volume for cap height |h| < c.
template <typename T> constexpr
T Geometric<T>::
ellipsoid_cap_vol(
    const T a,
    const T b,
    const T c,
    const T h
) noexcept
{    
    return pi<T> * a * b / (c*c) * h*h * (c - h / three<T>);
}


// Ellipsoidal cap base area for cap height |h| < c.
template <typename T> constexpr
T Geometric<T>::
ellipsoid_cap_base_area(
    const T a,
    const T b,
    const T c,
    const T h
) noexcept
{
    return pi<T> * a*b / (c*c) * h * (two<T> * c - h);
}


// Spheroidal surface area.
// Spheroid is given by r[0:2] = {a, b, c}, a = b,
// i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1
template <typename T> constexpr
T Geometric<T>::
spheroid_surf_area(
    const A3t& r,
    Msgr& msgr
)
{
    if (r[0] == r[1]) {
        const auto a2 = r[0]*r[0]; 
        const auto c2 = r[2]*r[2];
        
        if (r[0] < r[2]) {                    // prolate spheroid
            const auto e = std::sqrt(one<T> - a2 / c2);
            return twopi<T> * (a2 + std::asin(e) * r[0] * r[2] / e);
        }
        if (r[0] > r[2]) {                    // oblate spheroid
            const auto e = std::sqrt(one<T> - c2 / a2);
            return pi<T> * (two<T>*a2 +
                            std::log((one<T> + e) / (one<T> - e)) * c2 / e);
        }
        return four<T> * pi<T> * a2;        // sphere
    }
    throw SimpleException(
        "Error in spheroid_surf_area: spheroid r[0] == r[1] is required", &msgr);
}


// Determine spheroid axes symmetry from its dimensions r:
// returns (-1, -1, -1) if spheroid is a shpere,
// otherwise returns (i, j, k) where (i,j) are axes indexes of unequal
// dimensions and k is index of the pole axis.
template <typename T> inline
auto Geometric<T>::
spheroid_axes_symmetry(
    const A3t& r
) noexcept -> A3i
{
    if (all_sides_are_equal(r)) return A3i{-1};  // is a sphere

    if (r[0] == r[1]) return {0, 2, 2};    // pole is along 2
    if (r[0] == r[2]) return {0, 1, 1};    // pole is along 1
    if (r[1] == r[2]) return {0, 1, 0};    // pole is along 0

    XASSERT(false, " Geometric::spheroid_axes_symmetry failed");

    return {};
}


// Returns dimensions {a, b} of an ellipse x^2/a^2 + y^2/b^2 = 1
// resulting from the crossection of an ellipsoid
// x^2/e[0]^2 + y^2/e[1]^2 + z^2/e[2]^2 = 1
// with plane z = h.
template <typename T> constexpr
auto Geometric<T>::
ellipsoid_horizontal_crosection0(
    const A3t& e,
    const T h
) noexcept -> A2t
{
    const auto u = std::sqrt(e[2]*e[2] - h*h) / e[2];
    return { u * e[0],
             u * e[1] };
}


// Intersection of a line and an ellipsoid
// d: the line is given by a point v and a direction unit vector d
// e are the ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 semiaxes e = {a, b, c},
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipsoid(
    const A3t& v,
    const A3t& e,
    const A3t& d
) noexcept
{
    const auto u1 = d[0]*d[0] * e[1]*e[1] * e[2]*e[2] +
                    d[1]*d[1] * e[0]*e[0] * e[2]*e[2] +
                    d[2]*d[2] * e[0]*e[0] * e[1]*e[1];

    const auto u2 = ( d[0] * v[0] * e[1]*e[1] * e[2]*e[2] +
                      d[1] * v[1] * e[0]*e[0] * e[2]*e[2] +
                      d[2] * v[2] * e[0]*e[0] * e[1]*e[1] ) * two<T>;

    const auto u3 = v[0]*v[0] * e[1]*e[1] * e[2]*e[2] +
                    v[1]*v[1] * e[0]*e[0] * e[2]*e[2] +
                    v[2]*v[2] * e[0]*e[0] * e[1]*e[1] -
                    e[0]*e[0] * e[1]*e[1] * e[2]*e[2];

    const auto discr = u2 * u2 - four<T> * u1 * u3;
    
    if (discr >= zero<T>) {    // an intersection is possible

        // Putative intersection points:
        const auto t1 = (-u2 - std::sqrt(discr)) / (two<T> * u1);
        const auto t2 = (-u2 + std::sqrt(discr)) / (two<T> * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero<T> && t2 > zero<T>)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1:
        if (t1 > zero<T>) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero<T>) return t2;

        return huge<T>;     // no intersection
    }

    return huge<T>;    // no intersection
}


// Intersection of a line and a not rotated ellipse centered at the origin.
// The line is given by a point 'v' = (v0, v1) in plane and
// a direction vector 'd' = (d0, d1).
// The ellipse has dimensions 'ab' = (a, b): x^2 / a^2 + y^2 / b^2 = 1.
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipse(
    const A2t& v,
    const A2t& ab,
    const A2t& d
) noexcept
{
    const auto a2 = ab[0] * ab[0];
    const auto b2 = ab[1] * ab[1];

    const auto u1 = d[0] * d[0] * b2 +
                    d[1] * d[1] * a2;
    const auto u2 = two<T> * (d[0] * v[0] * b2 +
                              d[1] * v[1] * a2);
    const auto u3 = v[0] * v[0] * b2 +
                    v[1] * v[1] * a2 - a2 * b2;
    const auto discr = u2*u2 - four<T> * u1 * u3;
    
    if (discr >= zero<T>) {       // an intersection is possible
        const auto sd = std::sqrt(discr);
        // Putative intersection points:
        const auto t1 = (- u2 - sd) / (two<T> * u1);
        const auto t2 = (- u2 + sd) / (two<T> * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero<T> &&
            t2 > zero<T>)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1.
        if (t1 > zero<T>) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero<T>) return t2;
    }
    return huge<T>; // no intersection
}


// Intersection of a line and a rotated ellipse centered at the origin
// The line is given by a point 'v' = (v0, v1) in plane and a direction
// vector 'd' = (d0, d1). The ellipse is rotated counterclockwise through
// angle alpha about the origin, has semi-axes 'ab' = (a, b):
// (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1
// see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
template <typename T> constexpr
T Geometric<T>::
intersection_line_ellipse( const A2t& v,
                           const A2t& d,
                           const A2t& ab,
                           const T alpha ) noexcept
{
    const auto sia = std::sin(alpha);
    const auto coa = std::cos(alpha);
    const auto sia2 = sia*sia;
    const auto coa2 = coa*coa;

    const auto a2 = ab[0]*ab[0];
    const auto b2 = ab[1]*ab[1];

    // Coefficients of the quadratic form of the ellipse eq:
    // A*x^2 + B*x*y + C*y^2 = 1:
    const auto A = coa2 / a2 + sia2 / b2;
    const auto B = two<T> * coa * sia * (one<T> / a2 - one<T> / b2);
    const auto C = sia2 / a2 + coa2 / b2;

    // Coefficients of the quadratic equation of the line-ellipse intersection:
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

    if (discr >= zero<T>) {
        // There is an intersection possible:
        const auto sd = std::sqrt(discr);
        // Putative intersection points.
        const auto t1 = (- u2 - sd) / (two<T> * u1);
        const auto t2 = (- u2 + sd) / (two<T> * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2:
        if (t1 > zero<T> &&
            t2 > zero<T>)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1:
        if (t1 > zero<T>) return t1;

        // Only t2 is in the positive half-line: intersection is at t2:
        if (t2 > zero<T>) return t2;
    }
    return huge<T>;     // no intersection
}


/*
// Distance of a point p from a plane given by eq. dot(n,x) + o = 0
template <typename T> constexpr
T Geometric<T>::
distance_point_plane(
    const A3t& p,
    const A3t n,
    const A3t o
) noexcept
{
    return p.dot(n) + o;
}

// intersection of a line segment and and a plane
template <typename T> constexpr
bool Geometric<T>::
intersection_segment_plane(
    const A3t& p1,
    const A3t& p2,
    const A3t& n,
    const A3t& o,
    A3t& intersP
) noexcept
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


// Intersection of a line and and a plane.
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane(
    const A3t& pab,
    const A3t& p10,
    const A3t& p20,
    const A3t& pa0,
    const T s,
    Msgr& msgr
) noexcept
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


// Intersection of a line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by three points 'v1', 'v2', 'v3'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane(
    const A3t& p,
    const A3t& d,
    const A3t& v1,
    const A3t& v2,
    const A3t& v3
) noexcept
{
    const auto v13 = v1 - v3;
    const auto v23 = v2 - v3;

    const auto vn = A3t::crosspr(v13, v23);
    const auto n = vn.unitv();
    const auto nd = A3t::dotpr(n, d);

    // No intersection if the line is parallel to the plane:
    if (std::abs(nd) < EPS<T>)
        return -huge<T>;

    // Intersection point is at: p + t*d:
    const auto nw = -A3t::dotpr(n, p - v1);
    return nw / nd;
}


// Intersection of a line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by a points 'v', and a normal 'n'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template <typename T> constexpr
T Geometric<T>::
intersection_line_plane(
    const A3t& p,
    const A3t& d,
    const A3t& v,
    const A3t& n
) noexcept
{
    const auto nd = A3t::dotpr(n, d);

    // No intersection if the line is parallel to the plane:
    if (std::abs(nd) < EPS<T>)
        return -huge<T>;

    // Intersection point is at: p + t*d:
    return -A3t::dotpr(n, p - v) / nd;
}


// Intersection of a vertical line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by a points 'v', and a normal 'n'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template <typename T> constexpr
T Geometric<T>::
intersection_vertical_line_plane(
    const A3t& p,
    const int sign,
    const A3t& v,
    const A3t& n
) noexcept
{
    // No intersection if the line is parallel to the plane:
    if (std::abs(n[2]) < EPS<T>)
        return -huge<T>;

    // Intersection point is at: p + t*d:
    return -A3t::dotpr(n, p - v) / sign*n[2];
}


// Intersection of a line and and a cone.
template <typename T> constexpr
T Geometric<T>::
intersection_line_cone(
    const A3t& w,
    const A3t& q,
    const A3t& h,
    const A3t& m,
    const A3t& p,
    const A3t& d
) noexcept
{
    const auto u1 = h.dotpr(d*d) +
                    two<T>*(q[0]*q[1] * d[0]*d[1] +
                            q[0]*q[2] * d[0]*d[2] +
                            q[1]*q[2] * d[1]*d[2]);
                                     
    const auto u2 = m.dotpr(d) +
                    two<T>*(h.dotpr(d*p) +
                            q[0]*q[1] * (d[0]*p[1] + d[1]*p[0]) +
                            q[0]*q[2] * (d[0]*p[2] + d[2]*p[0]) +
                            q[1]*q[2] * (d[1]*p[2] + d[2]*p[1]));
                                                    
    const auto u3 = m.dotpr(p) +
                    h.dotpr(p*p + w*w) +
                    two<T>*(q[0]*q[1] * (p[0]*p[1] + w[0]*w[1]) +
                            q[0]*q[2] * (p[0]*p[2] + w[0]*w[2]) +
                            q[1]*q[2] * (p[1]*p[2] + w[1]*w[2]));

    auto discr = u2*u2 - four<T> * u1*u3;
    
    if (discr >= zero<T>) {
        // An intersection is possible.
        // Putative intersection points:
        const auto t1 = (- u2 - std::sqrt(discr)) / (two<T> * u1);
        const auto t2 = (- u2 + std::sqrt(discr)) / (two<T> * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero<T> &&
            t2 > zero<T>)
                 return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1.
        if (t1 > zero<T>) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero<T>) return t2;

        return -one<T>;     // No intersection.
    }

    return -one<T>;    // No intersection.
}


// Projection of a vector 'd' to z=0 plane.
template <typename T> constexpr
auto Geometric<T>::
proj2z0plane(
    const A3t& d
) noexcept -> A2t
{
    return {d(0,1).unitv(), zero<T>};
}


// Returns as 'rm' a rotation matrix for rotation over 'angle'
// around a general axis 'n'.
template <typename T> constexpr
void Geometric<T>::
rotmat(
    const A3t& n,
    T angle,
    T rm[3][3]
) noexcept
{
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    const auto c = - n*n + one<T>;
    
    rm[0][0] = n[0]*n[0] + c[0] * coa;
    rm[0][1] = n[0]*n[1] * (one<T> - coa) - n[2] * sia;
    rm[0][2] = n[0]*n[2] * (one<T> - coa) + n[1] * sia;

    rm[1][0] = n[0]*n[1] * (one<T> - coa) + n[2] * sia;
    rm[1][1] = n[1]*n[1] + c[1] * coa;
    rm[1][2] = n[1]*n[2] * (one<T> - coa) - n[0] * sia;

    rm[2][0] = n[0]*n[2] * (one<T> - coa) - n[1] * sia;
    rm[2][1] = n[1]*n[2] * (one<T> - coa) + n[0] * sia;
    rm[2][2] = n[2]*n[2] + c[2] * coa;
}


// Returns as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'x'.
template <typename T> constexpr
void Geometric<T>::
rotmatx(
    T angle,
    T rm[3][3]
) noexcept
{
    const A3t n {one<T>, zero<T>, zero<T>};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = n[0]*n[0] + (one<T> - n[0]*n[0]) * coa;
    rm[0][1] = zero<T>;
    rm[0][2] = zero<T>;

    rm[1][0] = zero<T>;
    rm[1][1] = zero<T>;
    rm[1][2] = - n[0] * sia;

    rm[2][0] = zero<T>;
    rm[2][1] = n[0] * sia;
    rm[2][2] = zero<T>;
}


// Returns as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'y'.
template <typename T> constexpr
void Geometric<T>::
rotmaty(
    T angle,
    T rm[3][3]
) noexcept
{
    const A3t n {zero<T>, one<T>, zero<T>};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = zero<T>;
    rm[0][1] = zero<T>;
    rm[0][2] = n[1] * sia;

    rm[1][0] = zero<T>;
    rm[1][1] = n[1]*n[1] + (one<T> - n[1]*n[1]) * coa;
    rm[1][2] = zero<T>;

    rm[2][0] = - n[1] * sia;
    rm[2][1] = zero<T>;
    rm[2][2] = zero<T>;
}


// Returns as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'z'.
template <typename T> constexpr
void Geometric<T>::
rotmatz(
    T angle,
    T rm[3][3]
) noexcept
{
    const A3t n {zero<T>, zero<T>, one<T>};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = zero<T>;
    rm[0][1] = - n[2] * sia;
    rm[0][2] = zero<T>;

    rm[1][0] = n[2] * sia;
    rm[1][1] = zero<T>;
    rm[1][2] = zero<T>;

    rm[2][0] = zero<T>;
    rm[2][1] = zero<T>;
    rm[2][2] = n[2]*n[2] + (one<T> - n[2]*n[2]) * coa;
}


// Unit normal at point p on surface of an axis-aligned ellipsoid
// (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
// with dimensions r = {a,b,c}.
template <typename T> constexpr
auto Geometric<T>::
unormal_on_ellipsoid(
    const A3t& r,
    const A3t& p
) noexcept -> A3t
{    
    const auto a2 = r[0] * r[0];
    const auto b2 = r[1] * r[1];
    const auto c2 = r[2] * r[2];
    const A3t n {b2 * c2 * p[0],
                   a2 * c2 * p[1],
                   a2 * b2 * p[2]};

    return n.unitv();
}


// Given the dimensions of a 3D body, determine if all_sides_are_equal.
template <typename T> constexpr
bool Geometric<T>::
all_sides_are_equal(
    const A3t& e
) noexcept
{
    return e[0] == e[1] &&
           e[1] == e[2];
}


// Given the dimensions of a 3D body, determine if two_sides_are_equal.
template <typename T> constexpr
bool Geometric<T>::
two_sides_are_equal(
    const A3t& e
) noexcept
{
    return e[0] == e[1] ||
           e[0] == e[2] ||
           e[1] == e[2];
}


template <typename T> constexpr
auto Geometric<T>::
sph2cart(
    const T ph,        // inclination
    const T th,        // azimuth
    const T rad
) -> A3t
{
    const auto coph = std::cos(ph);

    return A3t{
        coph * std::cos(th),
        coph * std::sin(th),
               std::sin(ph)
    } * rad;
}


// Conversion of spherical to cartesian coordinates.
// phi: inclination
// theta: azimuth
template <typename T> constexpr
auto Geometric<T>::
sphere2cart(
    const T r,
    const T theta,
    const T sinPhi,
    const T cosPhi
) noexcept -> A3t
{
    return { r * std::sin(theta) * sinPhi,
             r * std::cos(theta) * sinPhi,
             r * cosPhi };
}


// Conversion of polar to cartesian coordinates.
// phi: inclination
// theta: azimuth
template <typename T> constexpr
auto Geometric<T>::
polar2cart(
    const T r,
    const T theta
) noexcept -> A2t
{
    return { r * std::sin(theta),
             r * std::cos(theta) };
}


// Point on the line closest to the origin.
template <typename T> constexpr
auto Geometric<T>::
ptClosest2orgn(
    const A3t& p,
    const A3t& d
) noexcept -> A3t
{
    // https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    // A 3D line is defined with 6 Plücker coordinates L = ( d, p × d )
    // where d is the direction of the line, and p is any point along the line

    return A3t::crosspr(d, A3t::crosspr(d, p));
}

//----------------------------------------------------------------------------------------------------------------------

// point on an ellipsoid closest to line
template <typename T> constexpr
auto Geometric<T>::
ellipsoid_closest_point2Line(
    const A3t& ptc2o,
    const A3t& elps
) noexcept -> A3t
{
    // https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    // ptc2o is a point on a 3D line closest to the origin; 
    // center of the ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    // with dimensions elps = {a,b,c} is at the origin.

    const auto elps2 = elps * elps;
    const auto l = std::sqrt(A3t::dotpr(ptc2o*ptc2o, elps2));

    return (ptc2o * elps2) / l;
}


// Projection of vector v on a plane given by the normal n.
template <typename T> constexpr
auto Geometric<T>::
vector_proj2plane(
    const A3t& v,
    const A3t& n
) noexcept -> A3t
{
    return v - v.vecProjection(n);
}


// Squared distance between two line segments in 3D.
template <typename T> constexpr
T Geometric<T>::
squared_dist3D_Segment_to_Segment(
    const A3t& S10,
    const A3t& S11,
    const A3t& S20,
    const A3t& S21
) noexcept
{
    constexpr auto SMALL_NUM = static_cast<T>(0.00000001);   // anything that avoids division overflow
    const auto u = S11 - S10;
    const auto v = S21 - S20;
    const auto w = S10 - S20;
    const auto a = A3t::dotpr(u,u);            // always >= 0
    const auto b = A3t::dotpr(u,v);
    const auto c = A3t::dotpr(v,v);            // always >= 0
    const auto d = A3t::dotpr(u,w);
    const auto e = A3t::dotpr(v,w);
    const auto D = a*c - b*b;                    // always >= 0
    T sc, sN, sD = D;            // sc = sN / sD, default sD = D >= 0
    T tc, tN, tD = D;            // tc = tN / tD, default tD = D >= 0
    
    // Compute the line Config of the two closest points.
    if (D < SMALL_NUM) {
        // The lines are almost parallel.
        sN = zero<T>;         // force using point P0 on segment S1
        sD = one<T>;          // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {
        // Get the closest points on the infinite lines:
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < zero<T> ) {    // sc < 0 => the s=0 edge is visible
            sN = zero<T>;
            tN = e;
            tD = c;
        }
        else if (sN > sD ) {    // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
    if (tN < zero<T>) {
        // tc < 0 => the t=0 edge is visible:
        tN = zero<T>;
        // Recompute sc for this edge:
        if (-d < zero<T>) sN = zero<T>;
        else if (-d > a)  sN = sD;
        else {            sN = -d;
                          sD = a;
        }
    }
    else if (tN > tD) {
        // tc > 1  => the t=1 edge is visible:
        tN = tD;
        // Recompute sc for this edge:
        if ((-d + b) < zero<T>) sN = zero<T>;
        else if ((-d + b) > a)    sN = sD;
        else {                    sN = (-d +  b);
                                sD = a;
        }
    }
    // Finally do the division to get sc and tc:
    sc = (std::abs(sN) < SMALL_NUM ? zero<T> : sN / sD);
    tc = (std::abs(tN) < SMALL_NUM ? zero<T> : tN / tD);
    
    // Get the difference of the two closest points:
    const auto dP = w + (u * sc) - (v * tc);  // =  S1(sc) - S2(tc)
    
    return dP.dotpr();   // the closest distance
}


template <typename T> constexpr
auto Geometric<T>::
hexagonal_lattice(
    const A2t orig,
    const T step,
    const szt numLayers
) noexcept -> std::vector<A2t>
{
    std::vector<A2t> v;

    auto add_point = [&](const szt j, const int sign)
    {
        for (szt i=1; i<=2*numLayers+1-j; i++) {
            const auto a1 = half<T>*j + i - numLayers - one<T>;
            const auto a2 = sign * static_cast<T>(j) * std::sin(pi<T>/three<T>);
            v.emplace_back(orig + A2t {a1, a2} * step);
        }
    };

    add_point(0, 0);
    for (szt j=1; j<=numLayers; j++) {
        add_point(j, 1);
        add_point(j, -1);
    }

    return v;
}


template <typename T> constexpr
auto Geometric<T>::
numLayers_hexagonal_lattice (
    const szt numVertices
) noexcept -> szt
{
    constexpr const int N = 6;  // Number of sides in the hexagon.
    if (numVertices <= 1)
        return numVertices;

    szt numLayers {1};
    szt n {1};
    do n += N * numLayers++;
    while (n < numVertices);

    return numLayers;
}


template <typename T> constexpr
T Geometric<T>::
grad2rad(
    const T grad
) noexcept
{
    return grad * pi<T> / RAD2GRAD;
}


template <typename T> constexpr
T Geometric<T>::
rad2grad(
    const T rad
) noexcept
{
    return rad * RAD2GRAD / pi<T>;
}


template <typename T> constexpr
T Geometric<T>::
cos_two_segments(
    const A3t& p1,
    const A3t& p2,
    const A3t& p3,
    const A3t& p4
) noexcept
{
    const auto d1 = p2 - p1;
    const auto d2 = p4 - p3;

    return d1.dotpr(d2) / (d1.norm() * d2.norm());
}


template <typename T> constexpr
T Geometric<T>::
dotpr_two_segments(
    const A3t& p1,
    const A3t& p2,
    const A3t& p3,
    const A3t& p4
) noexcept
{
    const auto d1 = p2 - p1;
    const auto d2 = p4 - p3;

    return d1.dotpr(d2);
}


// Determines if a point is inside a 2D triangle.
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template <typename T> inline
bool Geometric<T>::
point_in_triangle(
    const A2t& p,
    const A2t& v1,
    const A2t& v2,
    const A2t& v3
) noexcept
{
    auto sign = [](
        const A2t& p1,
        const A2t& p2,
        const A2t& p3 ) noexcept
    {
        const auto d {(p1[0] - p3[0]) * (p2[1] - p3[1]) -
                      (p2[0] - p3[0]) * (p1[1] - p3[1])};
        return d;
    };

    const T d1 {sign(p, v1, v2)};
    const T d2 {sign(p, v2, v3)};
    const T d3 {sign(p, v3, v1)};

    return !((d1 < zero<T> || d2 < zero<T> || d3 < zero<T>) &&
             (d1 > zero<T> || d2 > zero<T> || d3 > zero<T>));
}


// Determines if a point is inside a 2D triangle.
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template <typename T> inline
bool Geometric<T>::
point_in_triangle(
    const T* p,
    const T* v1,
    const T* v2,
    const T* v3
) noexcept
{
    auto sign = [](
        const T* p1,
        const T* p2,
        const T* p3 ) noexcept
    {
        const auto d {(*p1 - *p3) * (*(p2+1) - *(p3+1)) -
                      (*p2 - *p3) * (*(p1+1) - *(p3+1))};
        return d;
    };

    const T d1 {sign(p, v1, v2)};
    const T d2 {sign(p, v2, v3)};
    const T d3 {sign(p, v3, v1)};

    return !((d1 < zero<T> || d2 < zero<T> || d3 < zero<T>) &&
             (d1 > zero<T> || d2 > zero<T> || d3 > zero<T>));
}

}  // namespace utils::common

#endif     // UTILS_COMMON_GEOMETRIC_FUNCTIONS
