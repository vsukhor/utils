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


#include <concepts>
#include <cmath>
#include <vector>

#include "../arrays/all.h"
#include "../msgr.h"
#include "../constants.h"
#include "misc.h"


namespace utils::common {

/**
 * \struct Geometric geometric_functions.h
 * \brief A loose collection of geometry-related static functions.
 * \tparam real Floating point type.
 */
template<std::floating_point real>
struct Geometric {

    using A2r = arrays::A2<real>;
    using A3r = arrays::A3<real>;
    using A3i = arrays::A3<int>;

    // Named real constants
    static constexpr auto zero = utils::zero<real>;
    static constexpr auto half = utils::half<real>;
    static constexpr auto one = utils::one<real>;
    static constexpr auto two = utils::two<real>;
    static constexpr auto three = utils::three<real>;
    static constexpr auto four = utils::four<real>;
    static constexpr auto pi = utils::pi<real>;
    static constexpr auto twopi = utils::twopi<real>;
    static constexpr auto halfpi = utils::halfpi<real>;
    static constexpr auto EPS = utils::EPS<real>;

    static constexpr auto RAD2GRAD = static_cast<real>(180);

    // Elliptic shapes +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Gets a point on an ellipse centered at zero.
     * \param alpha Angular coordinate.
     * \param ab    Dimensions of ellipse semi-major axes.
     */
    static A2r ellipse(
        real alpha,     
        const A2r& ab    
    ) noexcept;

    /**
     * \brief Point inside an ellipsoid or outside.
     * \details Finds out if point \p p is inside an ellipsoid given by
     * dimensions \p e of the ellipsoid semi-major axes.
     * \param p Point coordinates.
     * \param e Dimensions of the ellipsoid semi-major axes.
     */
    static constexpr bool is_inside_ellipsoid(
        const A3r& p,
        const A3r& e
    ) noexcept;
    
    /**
     * \brief Calculates area of an ellipse.
     * \param a Length of the semi-major axis.
     * \param b Length of the semi-major axis.
     */
    static constexpr real ellipse_area(
        real a, 
        real b
    ) noexcept;

    /**
     * \brief Calculates volume of an ellipsoid.
     * \param a Length of the semi-major axis.
     * \param b Length of the semi-major axis.
     * \param c Length of the semi-major axis.
     */
    static constexpr real ellipsoid_vol(
        real a, 
        real b, 
        real c
    ) noexcept;
     
    /**
     * \brief Calculates volume of an elliptic cylinder.
     * \param a Length of the semi-major axis.
     * \param b Length of the semi-major axis.
     * \param h Cylinder height.
     */
    static constexpr real elliptic_cylinder_vol(
        real a, 
        real b, 
        real h
    ) noexcept;
    
    /**
     * \brief Calculates volume of an ellipsoidal cap.
     * \param a Length of the semi-major axis.
     * \param b Length of the semi-major axis.
     * \param c Length of the semi-major axis.
     * \param h The cap height (|h| < c).
     */
    static constexpr real ellipsoid_cap_vol(
        real a, 
        real b, 
        real c, 
        real h
    ) noexcept;
     
    /**
     * \brief Calculates base area of an ellipsoidal cap.
     * \param a Length of the semi-major axis.
     * \param b Length of the semi-major axis.
     * \param c Length of the semi-major axis.
     * \param h The cap height (|h| < c).
     */
    static constexpr real ellipsoid_cap_base_area(
        real a, 
        real b,  
        real c, 
        real h
    ) noexcept;
    
    /**
     * \brief Calculates surface area of a spheroid.
     * \details Spheroid is given by \p r .
     * r  = {a, b, c}, a = b, i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1.
     * \param r Spheroid dimensions at semi-major axes.
     * \param msgr Output message processor.
     */
    static constexpr real spheroid_surf_area(
        const A3r& r,    
        Msgr& msgr       
    );

    /** 
     * \brief Unit normal on the surface of an axis-aligned ellipsoid.
     * \details Calculates unit normal vector at point \a p on surface
     * of an axis-aligned ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
     * with dimensions \p r = {a,b,c}.
     * \param r Dimensions of an ellipsoid.
     * \param p Point on ellipsoid surface.
     */
    static constexpr auto unormal_on_ellipsoid(
        const A3r& r,    
        const A3r& p     
        ) noexcept -> A3r;

    /** 
     * \brief Determines symmetry axes of a spheroid.
     * \details The spheroid should be axis-aligned.
     * \param r Dimensions of the spheroid semi-major axes.
     * \return (-1, -1, -1) if spheroid is a shpere,
     * otherwise (i, j, k) where (i,j) are axes indexes of unequal dimensions
     * and k is index of the pole axis.
     */
    static auto spheroid_axes_symmetry(
        const A3r& r
    ) noexcept -> A3i;

    /**
     * \brief Ellipse resulting from a horizontal cross-section of an ellipsoid.
     * \details Calculates dimensions {a, b} at semi-axes of an
     * ellipse x^2/a^2 + y^2/b^2 = 1 resulting from the horizontal z = h plane
     * cross-section of an ellipsoid x^2/e[0]^2 + y^2/e[1]^2 + z^2/e[2]^2 = 1
     * having dimensions \p e.
     * \param e Dimensions of the ellipsoid semi-axes.
     * \param h z-coordinate of the horizontal plane.
     * \return Semi-axes of the ellipsoid cross-section.
     */
    static constexpr auto ellipsoid_horizontal_crosection0(
        const A3r& e,   
        real h           
    ) noexcept -> A2r;


    // Line intersections ++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Intersection of a line and an ellipsoid.
     * \details Finds intersection of a line through a point \p v in the direction
     * unit vector \p d and an ellipsoid. The ellipsoid
     * x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 is given by its semi-axes e = {a, b, c}.
     * \param v Point on a line.
     * \param e Ellipsoid semi-axes.
     * \param d Direction of the line.
     */
    static constexpr real intersection_line_ellipsoid(
        const A3r& v,    
        const A3r& e,    
        const A3r& d     
    ) noexcept;

    /**
     * \brief Intersection of a line and an ellipse.
     * Find intersection of a line and a not rotated ellipse centered at
     * the origin. The line is given by a point \p v in plane and a
     * direction vector \p d = (d0, d1) The ellipse has dimensions
     * \p e = (a, b): x^2 / a^2 + y^2 / b^2 = 1.
     * \param v  Point on the line.
     * \param ab Dimensions of the ellipse semi-major axes.
     * \param d  Line direction vector.
     */
    static constexpr real intersection_line_ellipse(
        const A2r& v,     
        const A2r& ab,    
        const A2r& d      
    ) noexcept;

    /**
     * \brief Intersection of a line and a rotated ellipse centered at the origin.
     * \details The line is given by a point \p v = (v0, v1) in plane and a
     * direction vector \p d = (d0, d1). The ellipse is rotated counterclockwise
     * over angle alpha about the origin, has semi-axes \p ab = (a, b)
     * (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1.
     * \see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
     * \param v  Point on the line.
     * \param d  Line direction vector.
     * \param ab Dimensions of the ellipse semi-major axes.
     */
    static constexpr real intersection_line_ellipse(
        const A2r& v,    
        const A2r& d,    
        const A2r& ab,   
        real alpha
    ) noexcept;

    // Distance of a point p from a plane given by eq. dot(n,x) + o = 0.
    // static constexpr real distance_point_plane(
    //      const A3r& p, const A3r n, const A3r o ) noexcept;

    // Intersection of a line segment and a plane.
    // static constexpr bool intersection_segment_plane(
    //      const A3r& p1, const A3r& p2, const A3r& n,
    //      const A3r& o, A3r& intersP ) noexcept;
    
    /**
     * \brief Finds intersection of a line and a plane.
     */
    static constexpr real intersection_line_plane(
        const A3r& pab,
        const A3r& p10,
        const A3r& p20,
        const A3r& pa0,
        real s,
        Msgr &msgr
    ) noexcept;

    /**
     * \brief Finds intersection of a line and a plane.
     * \details The line is defined by a point \p p0 and direction vector \p d .
     * The plane is defined by three points \p v1  \p v2  \p v3 .
     * \param p0 Point on the line.
     * \param d  Line direction vector.
     * \param v1 Point on the plane.
     * \param v2 Point on the plane.
     * \param v3 Point on the plane.
     * \return Distance in direction \p d from \p p0 to the intersection point.
     */
    static constexpr real intersection_line_plane(
        const A3r& p0,      
        const A3r& d,        
        const A3r& v1,       
        const A3r& v2,       
        const A3r& v3        
    ) noexcept;

    /**
     * \brief Finds intersection of a line and a plane.
     * \details The line is defined by a point \p p0 and direction vector \p d .
     * The plane is defined by a point \p v and a normal \p n .
     * \param p0 Point on the line.
     * \param d  Line direction vector.
     * \param v  Point on the plane.
     * \param n  Plane unit normal vector.
     * \return Distance in direction \p d from \p p0 to the intersection point.
     */
    static constexpr real intersection_line_plane(
        const A3r& p0,       
        const A3r& d,        
        const A3r& v,        
        const A3r& n       
    ) noexcept;

    /**
     * \brief Finds intersection of a vertical line and a plane.
     * \details The line is defined by a point \p p0 .
     * The plane is defined by a point \p v, and a unit normal vector \p n .
     * \param p0   Point on the line.
     * \param sign Directionality (-1, 1) of the result relative to \p d .
     * \param v    Point on the plane.
     * \param n    Plane unit normal vector.
     * \return Distance between \p p0 and the intersection point; parallel or
     * antiparallel to the line -- depending on \p sign .
     */
    static constexpr real intersection_vertical_line_plane(
        const A3r& p0,   
        int sign,        
        const A3r& v,   
        const A3r& n    
    ) noexcept;

    /**
     * \brief Finds intersection of a line and a cone.
     */
    static constexpr real intersection_line_cone(
        const A3r& w,
        const A3r& q,
        const A3r& h,
        const A3r& m,
        const A3r& p,
        const A3r& d
    ) noexcept;


    // Rotations +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Rotation matrix around a general axis \p n .
     * \details Sets rotatios over \p angle around a general axis \p n .
     * \param[in] n     Direction of rotation axis.
     * \param[in] angle Rotation angle.
     * \param[out] rm   Rotation matrix.
     */
    static constexpr void rotmat(
        const A3r& n,
        real angle,
        real rm[3][3]
    ) noexcept;

    /**
     * \brief Rotation matrix around an axis parallel to 'x'.
     * \details Sets rotations over \p angle around an axis parallel to 'x'.
     * \param[in] angle Rotation angle.
     * \param[out] rm   Rotation matrix.
     */
    static constexpr void rotmatx(
        real angle,
        real rm[3][3]
    ) noexcept;

    /** 
     * \brief Rotation matrix around an axis parallel to 'y'.
     * \details Sets rotations over \p angle around an axis parallel to 'y'.
     * \param[in] angle Rotation angle.
     * \param[out] rm   Rotation matrix.
     */
    static constexpr void rotmaty(
        real angle,
        real rm[3][3]
    ) noexcept;

    /**
     * \brief Rotation matrix around an axis parallel to 'z'.
     * \details Sets rotations over \p angle around an axis parallel to 'z'.
     * \param[in] angle Rotation angle.
     * \param[out] rm   Rotation matrix.
     */
    static constexpr void rotmatz(
        real angle,
        real rm[3][3]
    ) noexcept;


    // Comparisons +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Determines if all sides of a 3D body are equal.
     * \details Given the dimensions \p e of a 3D body, determines if all its 
     * sides are equal.
     * \param e Body dimensions.
     */
    static constexpr bool all_sides_are_equal(const A3r& e) noexcept;
    
    /**
     * \brief Determines if two sides of a 3D body are equal.
     * \details Given the dimensions \p e of a 3D body, determines if two of its 
     * sides are equal.
     * \param e Body dimensions.
     */
    static constexpr bool two_sides_are_equal(const A3r& e) noexcept;
    

    // Some conversions ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Converts spherical coordinates to cartesian coordinates.
     * \param ph  Inclination.
     * \param th  Azimuth.
     * \param rad Radius.
     * \return (x, y, z).
     */
    static constexpr auto sph2cart(
        real ph,            
        real th,           
        real rad=one     
    ) -> A3r;

    /**
     * \brief Converts spherical coordinates to cartesian coordinates.
     * \param r Radius.
     * \return (x, y, z).
     */
    static constexpr auto sphere2cart(
        real r,        
        real theta,
        real sinPhi,
        real cosPhi
    ) noexcept  -> A3r;

    /**
     * \brief Converts polar coordinates to cartesian coordinates.
     * \param r Radius.
     * \param theta Angle.
     * \return (x, y).
     */
    static constexpr auto polar2cart(
        real r,        
        real theta     
    ) noexcept -> A2r;

    /**
     * \brief Converts Grad to Rad.
     * \param grad Angle given in Grads.
     */
    static constexpr real grad2rad(real grad) noexcept;

    /** 
     * \brief Converts Rad to Grad.
     * \param rad Angle given in Rads.
     */
    static constexpr real rad2grad(real rad) noexcept;


    // Closest points ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Finds point on a line closest to the origin.
     * Find point closest to the origin on the line that passes
     * through point \p p in the direction \p d . The 3D line is
     * defined with 6 Plücker coordinates L = (d, p × d),
     * where \p d is the direction of the line, and \p p is any point along the line.
     * \see https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
     * \param p Arbitrary point along the line.
     * \param d Direction of the line.
     * \return Point on the line closest to the origin.
     */
    static constexpr auto ptClosest2orgn(
        const A3r& p,    
        const A3r& d    
    ) noexcept -> A3r;
    
    /**
     * \brief Finds point on an ellipsoid closest to a line.
     * \details Finds point closest to line on an ellipsoid given by its
     * semi-major axes \p e . Center of the ellipsoid
     *      (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
     * with dimensions e = {a,b,c} is at the origin.
     * https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
     * \param ptc2o Point on a 3D line closest to the origin.
     * \param e     Ellipsoid given by its semi-major axes.
     * \return Point on an ellipsoid closest to line.
     */
    static constexpr auto ellipsoid_closest_point2Line(
        const A3r& ptc2o,    
        const A3r& e         
    ) noexcept  -> A3r;
    

    // Projections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Finds projection of a vector \b d to z=0 plane.
     */
    static constexpr auto proj2z0plane( const A3r& d ) noexcept -> A2r;

    /**
     * \brief Projection of vector \b v on a plane defined by the normal \b n .
     */
    static constexpr auto vector_proj2plane(
        const A3r& v,
        const A3r& n
    ) noexcept -> A3r;
    

    // Hexagonal lattice +++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Find coordinates of a hexagonal lattice.
     * \details Find coordinates of a hexagonal lattice centered at \p orig
     * with \p step and number of layers \p numLayers .
     */
    static constexpr auto hexagonal_lattice(
        A2r orig,
        real step,
        szt numLayers
    ) noexcept -> std::vector<A2r>;

    /**
     * \brief Finds number of layers in a hexagonal lattice.
     * \details The number of layers is determined from the number of vertices
     * \p numVertices .
     */
    static constexpr auto numLayers_hexagonal_lattice(
        szt numVertices
    ) noexcept -> szt;


    // Two segments ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Calculates cosine of an angle between two 3D segments.
     * \details The segments are defined by their end points.
     */
    static constexpr real cos_two_segments(
        const A3r& p1,
        const A3r& p2,
        const A3r& p3,
        const A3r& p4
    ) noexcept;

    /**
     * \brief Calculates dot product of two 3D vectors given by line segments.
     * \details The segments are defined by their end points.
     */
    static constexpr real dotpr_two_segments(
        const A3r& p1,
        const A3r& p2,
        const A3r& p3,
        const A3r& p4
    ) noexcept;

    /** 
     * \brief Calculates squared distance between two line segments in 3D.
     * \details The segments are defined by pairs of points
     * [ \p S10, \p S11 ] and [ \p S20, \p S21 ].
     */
    static constexpr real squared_dist3D_Segment_to_Segment(
        const A3r& S10, const A3r& S11,
        const A3r& S20, const A3r& S21
    ) noexcept;

    // Points inside triangle ++++++++++++++++++++++++++++++++++++++++++++++++++

    /**
     * \brief Determines if a point is inside a 2D triangle.
     * \details Find out if point \p pt is inside triangle given by its
     * vertexes \p v1, \p v2 and \p v3.
     * https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
     */
    static bool point_in_triangle(
        const A2r& pt,
        const A2r& v1,
        const A2r& v2,
        const A2r& v3
    ) noexcept;

    /**
     * \brief Determines if a point is inside a 2D triangle.
     * \details Finds out if point \p p is inside triangle given by its
     * vertexes \p v1, \p v2 and \p v3.
     * https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
     */
    static bool point_in_triangle (
        const real* p,
        const real* v1,
        const real* v2,
        const real* v3
    ) noexcept;
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


template<std::floating_point real>
auto Geometric<real>::
ellipse(
    const real alpha,
    const A2r& ab
) noexcept -> A2r
{
    return { ab[0]*std::cos(alpha),
             ab[1]*std::sin(alpha) };
}


// ellipse area
template<std::floating_point real> constexpr
real Geometric<real>::
ellipse_area(
    const real a,
    const real b
) noexcept
{
    return pi*a*b;
}


// ellipsoidal volume
template<std::floating_point real> constexpr
real Geometric<real>::
ellipsoid_vol(
    const real a,
    const real b,
    const real c
) noexcept
{
    return four/three * pi*a*b*c;
}


template<std::floating_point real> constexpr
bool Geometric<real>::
is_inside_ellipsoid(
    const A3r& p,
    const A3r& e
) noexcept
{
    return p[0]*p[0]/e[0]/e[0] +
           p[1]*p[1]/e[1]/e[1] +
           p[2]*p[2]/e[2]/e[2] < one - EPS;
}


// Elliptic cylinder volume.
template<std::floating_point real> constexpr
real Geometric<real>::
elliptic_cylinder_vol(
    const real a,
    const real b,
    const real h
) noexcept
{
    return h * ellipseArea(a, b);
}


// Ellipsoidal cap volume for cap height |h| < c.
template<std::floating_point real> constexpr
real Geometric<real>::
ellipsoid_cap_vol(
    const real a,
    const real b,
    const real c,
    const real h
) noexcept
{    
    return pi * a * b / (c*c) * h*h * (c - h / three);
}


// Ellipsoidal cap base area for cap height |h| < c.
template<std::floating_point real> constexpr
real Geometric<real>::
ellipsoid_cap_base_area(
    const real a,
    const real b,
    const real c,
    const real h
) noexcept
{
    return pi * a*b / (c*c) * h * (two * c - h);
}


// Spheroidal surface area.
// Spheroid is given by r[0:2] = {a, b, c}, a = b,
// i.e. (x^2+y^2)/a^2 + z^2/c^2 = 1
template<std::floating_point real> constexpr
real Geometric<real>::
spheroid_surf_area(
    const A3r& r,
    Msgr& msgr
)
{
    if (r[0] == r[1]) {
        const auto a2 = r[0]*r[0]; 
        const auto c2 = r[2]*r[2];
        
        if (r[0] < r[2]) {            // prolate spheroid
            const auto e = std::sqrt(one - a2 / c2);
            return twopi * (a2 + std::asin(e) * r[0] * r[2] / e);
        }
        if (r[0] > r[2]) {            // oblate spheroid
            const auto e = std::sqrt(one - c2 / a2);
            return pi * (two*a2 + std::log((one + e) / (one - e)) * c2 / e);
        }
        return four * pi * a2;        // sphere
    }
    throw common::Exception(
        "Error in spheroid_surf_area: spheroid r[0] == r[1] is required", &msgr);
}


// Determine spheroid axes symmetry from its dimensions r:
// returns (-1, -1, -1) if spheroid is a shpere,
// otherwise returns (i, j, k) where (i,j) are axes indexes of unequal
// dimensions and k is index of the pole axis.
template<std::floating_point real>
auto Geometric<real>::
spheroid_axes_symmetry(
    const A3r& r
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
template<std::floating_point real> constexpr
auto Geometric<real>::
ellipsoid_horizontal_crosection0(
    const A3r& e,
    const real h
) noexcept -> A2r
{
    const auto u = std::sqrt(e[2]*e[2] - h*h) / e[2];
    return { u * e[0],
             u * e[1] };
}


// Intersection of a line and an ellipsoid
// d: the line is given by a point v and a direction unit vector d
// e are the ellipsoid x^2/a^2 + y^2/b^2 + z^2/c^2 = 1 semiaxes e = {a, b, c},
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_ellipsoid(
    const A3r& v,
    const A3r& e,
    const A3r& d
) noexcept
{
    const auto u1 = d[0]*d[0] * e[1]*e[1] * e[2]*e[2] +
                    d[1]*d[1] * e[0]*e[0] * e[2]*e[2] +
                    d[2]*d[2] * e[0]*e[0] * e[1]*e[1];

    const auto u2 = (d[0] * v[0] * e[1]*e[1] * e[2]*e[2] +
                     d[1] * v[1] * e[0]*e[0] * e[2]*e[2] +
                     d[2] * v[2] * e[0]*e[0] * e[1]*e[1]) * two;

    const auto u3 = v[0]*v[0] * e[1]*e[1] * e[2]*e[2] +
                    v[1]*v[1] * e[0]*e[0] * e[2]*e[2] +
                    v[2]*v[2] * e[0]*e[0] * e[1]*e[1] -
                    e[0]*e[0] * e[1]*e[1] * e[2]*e[2];

    const auto discr = u2 * u2 - four * u1 * u3;
    
    if (discr >= zero) {    // an intersection is possible

        // Putative intersection points:
        const auto t1 = (-u2 - std::sqrt(discr)) / (two * u1);
        const auto t2 = (-u2 + std::sqrt(discr)) / (two * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero && t2 > zero)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1:
        if (t1 > zero) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero) return t2;

        return huge<real>;     // no intersection
    }

    return huge<real>;    // no intersection
}


// Intersection of a line and a not rotated ellipse centered at the origin.
// The line is given by a point 'v' = (v0, v1) in plane and
// a direction vector 'd' = (d0, d1).
// The ellipse has dimensions 'ab' = (a, b): x^2 / a^2 + y^2 / b^2 = 1.
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_ellipse(
    const A2r& v,
    const A2r& ab,
    const A2r& d
) noexcept
{
    const auto a2 = ab[0] * ab[0];
    const auto b2 = ab[1] * ab[1];

    const auto u1 = d[0] * d[0] * b2 +
                    d[1] * d[1] * a2;
    const auto u2 = two * (d[0] * v[0] * b2 +
                           d[1] * v[1] * a2);
    const auto u3 = v[0] * v[0] * b2 +
                    v[1] * v[1] * a2 - a2 * b2;
    const auto discr = u2*u2 - four * u1 * u3;
    
    if (discr >= zero) {       // an intersection is possible
        const auto sd = std::sqrt(discr);
        // Putative intersection points:
        const auto t1 = (- u2 - sd) / (two * u1);
        const auto t2 = (- u2 + sd) / (two * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero &&
            t2 > zero)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1.
        if (t1 > zero) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero) return t2;
    }
    return huge<real>; // no intersection
}


// Intersection of a line and a rotated ellipse centered at the origin
// The line is given by a point 'v' = (v0, v1) in plane and a direction
// vector 'd' = (d0, d1). The ellipse is rotated counterclockwise through
// angle alpha about the origin, has semi-axes 'ab' = (a, b):
// (x*cos(alpha) + ysin(alpha))^2 / a^2 + (x*sin(alpha) - y*cos(alpha))^2 / b^2 = 1
// see https://www.maa.org/external_archive/joma/Volume8/Kalman/General.html
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_ellipse( const A2r& v,
                           const A2r& d,
                           const A2r& ab,
                           const real alpha ) noexcept
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
    const auto B = two * coa * sia * (one / a2 - one / b2);
    const auto C = sia2 / a2 + coa2 / b2;

    // Coefficients of the quadratic equation of the line-ellipse intersection:
    const auto u1 = d[0]*d[0] * A +
                    d[0]*d[1] * B +
                    d[1]*d[1] * C;
    const auto u2 = two * (A * v[0] * d[0] +
                           C * v[1] * d[1]) +
                    B * (v[0] * d[1] +
                         v[1] * d[0]);
    const auto u3 = A * v[0]*v[0] +
                    B * v[0]*v[1] +
                    C * v[1]*v[1] - one;
    const auto discr = u2*u2 - four * u1 * u3;

    if (discr >= zero) {
        // There is an intersection possible:
        const auto sd = std::sqrt(discr);
        // Putative intersection points.
        const auto t1 = (- u2 - sd) / (two * u1);
        const auto t2 = (- u2 + sd) / (two * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2:
        if (t1 > zero &&
            t2 > zero)
            return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1:
        if (t1 > zero) return t1;

        // Only t2 is in the positive half-line: intersection is at t2:
        if (t2 > zero) return t2;
    }
    return huge<real>;     // no intersection
}


/*
// Distance of a point p from a plane given by eq. dot(n,x) + o = 0
template<std::floating_point real> constexpr
real Geometric<real>::
distance_point_plane(
    const A3r& p,
    const A3r n,
    const A3r o
) noexcept
{
    return p.dot(n) + o;
}

// intersection of a line segment and and a plane
template<std::floating_point real> constexpr
bool Geometric<real>::
intersection_segment_plane(
    const A3r& p1,
    const A3r& p2,
    const A3r& n,
    const A3r& o,
    A3r& intersP
) noexcept
{
    const auto d1 = distFromPlane(p1, n, o);
    const auto d2 = distFromPlane(p2, n, o);

    const bool bP1OnPlane = (std::abs(d1) < EPS);
    const bool bP2OnPlane = (std::abs(d2) < EPS);

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
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_plane(
    const A3r& pab,
    const A3r& p10,
    const A3r& p20,
    const A3r& pa0,
    const real s,
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

    return (p20[0] * (pa0[2] - p10[2]*lamda) -
            p20[2] * (pa0[0] - p10[0]*lamda)) / den2 * s;
}


// Intersection of a line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by three points 'v1', 'v2', 'v3'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_plane(
    const A3r& p,
    const A3r& d,
    const A3r& v1,
    const A3r& v2,
    const A3r& v3
) noexcept
{
    const auto v13 = v1 - v3;
    const auto v23 = v2 - v3;

    const auto vn = A3r::crosspr(v13, v23);
    const auto n = vn.unitv();
    const auto nd = A3r::dotpr(n, d);

    // No intersection if the line is parallel to the plane:
    if (std::abs(nd) < EPS)
        return -huge<real>;

    // Intersection point is at: p + t*d:
    const auto nw = -A3r::dotpr(n, p - v1);
    return nw / nd;
}


// Intersection of a line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by a points 'v', and a normal 'n'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_plane(
    const A3r& p,
    const A3r& d,
    const A3r& v,
    const A3r& n
) noexcept
{
    const auto nd = A3r::dotpr(n, d);

    // No intersection if the line is parallel to the plane:
    if (std::abs(nd) < EPS)
        return -huge<real>;

    // Intersection point is at: p + t*d:
    return -A3r::dotpr(n, p - v) / nd;
}


// Intersection of a vertical line and and a plane.
// The line is defined by a point 'p' and direction vector 'd'.
// The plane is defined by a points 'v', and a normal 'n'.
// Returns distance in direction 'd' from 'p0' to the intersection point.
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_vertical_line_plane(
    const A3r& p,
    const int sign,
    const A3r& v,
    const A3r& n
) noexcept
{
    // No intersection if the line is parallel to the plane:
    if (std::abs(n[2]) < EPS)
        return -huge<real>;

    // Intersection point is at: p + t*d:
    return -A3r::dotpr(n, p - v) / sign*n[2];
}


// Intersection of a line and and a cone.
template<std::floating_point real> constexpr
real Geometric<real>::
intersection_line_cone(
    const A3r& w,
    const A3r& q,
    const A3r& h,
    const A3r& m,
    const A3r& p,
    const A3r& d
) noexcept
{
    const auto u1 = h.dotpr(d*d) +
                    two*(q[0]*q[1] * d[0]*d[1] +
                            q[0]*q[2] * d[0]*d[2] +
                            q[1]*q[2] * d[1]*d[2]);
                                     
    const auto u2 = m.dotpr(d) +
                    two*(h.dotpr(d*p) +
                            q[0]*q[1] * (d[0]*p[1] + d[1]*p[0]) +
                            q[0]*q[2] * (d[0]*p[2] + d[2]*p[0]) +
                            q[1]*q[2] * (d[1]*p[2] + d[2]*p[1]));
                                                    
    const auto u3 = m.dotpr(p) +
                    h.dotpr(p*p + w*w) +
                    two*(q[0]*q[1] * (p[0]*p[1] + w[0]*w[1]) +
                            q[0]*q[2] * (p[0]*p[2] + w[0]*w[2]) +
                            q[1]*q[2] * (p[1]*p[2] + w[1]*w[2]));

    auto discr = u2*u2 - four * u1*u3;
    
    if (discr >= zero) {
        // An intersection is possible.
        // Putative intersection points:
        const auto t1 = (- u2 - std::sqrt(discr)) / (two * u1);
        const auto t2 = (- u2 + std::sqrt(discr)) / (two * u1);

        // Both t1 and t2 are in the positive half-line:
        // intersection is at the closest of t1, t2.
        if (t1 > zero &&
            t2 > zero)
                 return (t1 <= t2) ? t1 : t2;

        // Only t1 is in the positive half-line: intersection is at t1.
        if (t1 > zero) return t1;

        // Only t2 is in the positive half-line: intersection is at t2.
        if (t2 > zero) return t2;

        return -one;    // No intersection.
    }

    return -one;    // No intersection.
}


// Projection of a vector 'd' to z=0 plane.
template<std::floating_point real> constexpr
auto Geometric<real>::
proj2z0plane(
    const A3r& d
) noexcept -> A2r
{
    return {d(0,1).unitv(), zero};
}


// Return as 'rm' a rotation matrix for rotation over 'angle'
// around a general axis 'n'.
template<std::floating_point real> constexpr
void Geometric<real>::
rotmat(
    const A3r& n,
    real angle,
    real rm[3][3]
) noexcept
{
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    const auto c = - n*n + one;
    
    rm[0][0] = n[0]*n[0] + c[0] * coa;
    rm[0][1] = n[0]*n[1] * (one - coa) - n[2] * sia;
    rm[0][2] = n[0]*n[2] * (one - coa) + n[1] * sia;

    rm[1][0] = n[0]*n[1] * (one - coa) + n[2] * sia;
    rm[1][1] = n[1]*n[1] + c[1] * coa;
    rm[1][2] = n[1]*n[2] * (one - coa) - n[0] * sia;

    rm[2][0] = n[0]*n[2] * (one - coa) - n[1] * sia;
    rm[2][1] = n[1]*n[2] * (one - coa) + n[0] * sia;
    rm[2][2] = n[2]*n[2] + c[2] * coa;
}


// Return as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'x'.
template<std::floating_point real> constexpr
void Geometric<real>::
rotmatx(
    real angle,
    real rm[3][3]
) noexcept
{
    const A3r n {one, zero, zero};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = n[0]*n[0] + (one - n[0]*n[0]) * coa;
    rm[0][1] = zero;
    rm[0][2] = zero;

    rm[1][0] = zero;
    rm[1][1] = zero;
    rm[1][2] = - n[0] * sia;

    rm[2][0] = zero;
    rm[2][1] = n[0] * sia;
    rm[2][2] = zero;
}


// Return as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'y'.
template<std::floating_point real> constexpr
void Geometric<real>::
rotmaty(
    real angle,
    real rm[3][3]
) noexcept
{
    const A3r n {zero, one, zero};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = zero;
    rm[0][1] = zero;
    rm[0][2] = n[1] * sia;

    rm[1][0] = zero;
    rm[1][1] = n[1]*n[1] + (one - n[1]*n[1]) * coa;
    rm[1][2] = zero;

    rm[2][0] = - n[1] * sia;
    rm[2][1] = zero;
    rm[2][2] = zero;
}


// Return as 'rm' a rotation matrix for rotation over 'angle'
// around an axis parallel to 'z'.
template<std::floating_point real> constexpr
void Geometric<real>::
rotmatz(
    real angle,
    real rm[3][3]
) noexcept
{
    const A3r n {zero, zero, one};
    
    const auto sia = std::sin(angle);
    const auto coa = std::cos(angle);
    
    rm[0][0] = zero;
    rm[0][1] = - n[2] * sia;
    rm[0][2] = zero;

    rm[1][0] = n[2] * sia;
    rm[1][1] = zero;
    rm[1][2] = zero;

    rm[2][0] = zero;
    rm[2][1] = zero;
    rm[2][2] = n[2]*n[2] + (one - n[2]*n[2]) * coa;
}


// Unit normal at point p on surface of an axis-aligned ellipsoid
// (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
// with dimensions r = {a,b,c}.
template<std::floating_point real> constexpr
auto Geometric<real>::
unormal_on_ellipsoid(
    const A3r& r,
    const A3r& p
) noexcept -> A3r
{    
    const auto a2 = r[0] * r[0];
    const auto b2 = r[1] * r[1];
    const auto c2 = r[2] * r[2];

    const A3r n {b2 * c2 * p[0],
                 a2 * c2 * p[1],
                 a2 * b2 * p[2]};

    return n.unitv();
}


// Given the dimensions of a 3D body, determine if all_sides_are_equal.
template<std::floating_point real> constexpr
bool Geometric<real>::
all_sides_are_equal(
    const A3r& e
) noexcept
{
    return e[0] == e[1] &&
           e[1] == e[2];
}


// Given the dimensions of a 3D body, determine if two_sides_are_equal.
template<std::floating_point real> constexpr
bool Geometric<real>::
two_sides_are_equal(
    const A3r& e
) noexcept
{
    return e[0] == e[1] ||
           e[0] == e[2] ||
           e[1] == e[2];
}


template<std::floating_point real> constexpr
auto Geometric<real>::
sph2cart(
    const real ph,        // inclination
    const real th,        // azimuth
    const real rad
) -> A3r
{
    const auto coph = std::cos(ph);

    return A3r{
        coph * std::cos(th),
        coph * std::sin(th),
               std::sin(ph)
    } * rad;
}


// Conversion of spherical to cartesian coordinates.
// phi: inclination
// theta: azimuth
template<std::floating_point real> constexpr
auto Geometric<real>::
sphere2cart(
    const real r,
    const real theta,
    const real sinPhi,
    const real cosPhi
) noexcept -> A3r
{
    return { r * std::sin(theta) * sinPhi,
             r * std::cos(theta) * sinPhi,
             r * cosPhi };
}


// Conversion of polar to cartesian coordinates.
// phi: inclination
// theta: azimuth
template<std::floating_point real> constexpr
auto Geometric<real>::
polar2cart(
    const real r,
    const real theta
) noexcept -> A2r
{
    return { r * std::sin(theta),
             r * std::cos(theta) };
}


// Point on the line closest to the origin.
template<std::floating_point real> constexpr
auto Geometric<real>::
ptClosest2orgn(
    const A3r& p,
    const A3r& d
) noexcept -> A3r
{
    // https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    // A 3D line is defined with 6 Plücker coordinates L = ( d, p × d )
    // where d is the direction of the line, and p is any point along the line

    return A3r::crosspr(d, A3r::crosspr(d, p));
}

//----------------------------------------------------------------------------------------------------------------------

// point on an ellipsoid closest to line
template<std::floating_point real> constexpr
auto Geometric<real>::
ellipsoid_closest_point2Line(
    const A3r& ptc2o,
    const A3r& elps
) noexcept -> A3r
{
    // https://math.stackexchange.com/questions/895385/point-on-an-ellipsoid-closest-to-line?noredirect=1&lq=1
    // ptc2o is a point on a 3D line closest to the origin; 
    // center of the ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    // with dimensions elps = {a,b,c} is at the origin.

    const auto elps2 = elps * elps;
    const auto l = std::sqrt(A3r::dotpr(ptc2o*ptc2o, elps2));

    return (ptc2o * elps2) / l;
}


// Projection of vector v on a plane given by the normal n.
template<std::floating_point real> constexpr
auto Geometric<real>::
vector_proj2plane(
    const A3r& v,
    const A3r& n
) noexcept -> A3r
{
    return v - v.vecProjection(n);
}


// Squared distance between two line segments in 3D.
template<std::floating_point real> constexpr
real Geometric<real>::
squared_dist3D_Segment_to_Segment(
    const A3r& S10,
    const A3r& S11,
    const A3r& S20,
    const A3r& S21
) noexcept
{
    constexpr auto SMALL_NUM = static_cast<real>(0.00000001);   // anything that avoids division overflow
    const auto u = S11 - S10;
    const auto v = S21 - S20;
    const auto w = S10 - S20;
    const auto a = A3r::dotpr(u,u);            // always >= 0
    const auto b = A3r::dotpr(u,v);
    const auto c = A3r::dotpr(v,v);            // always >= 0
    const auto d = A3r::dotpr(u,w);
    const auto e = A3r::dotpr(v,w);
    const auto D = a*c - b*b;                    // always >= 0
    real sc, sN, sD = D;   // sc = sN / sD, default sD = D >= 0
    real tc, tN, tD = D;   // tc = tN / tD, default tD = D >= 0
    
    // Compute the line Config of the two closest points.
    if (D < SMALL_NUM) {
        // The lines are almost parallel.
        sN = zero;         // force using point P0 on segment S1
        sD = one;          // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {
        // Get the closest points on the infinite lines:
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < zero) {    // sc < 0 => the s=0 edge is visible
            sN = zero;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }
    if (tN < zero) {
        // tc < 0 => the t=0 edge is visible:
        tN = zero;
        // Recompute sc for this edge:
        if (-d < zero)   sN = zero;
        else if (-d > a) sN = sD;
        else {           sN = -d;
                         sD = a;
        }
    }
    else if (tN > tD) {
        // tc > 1  => the t=1 edge is visible:
        tN = tD;
        // Recompute sc for this edge:
        if      ((-d + b) < zero) sN = zero;
        else if ((-d + b) > a)    sN = sD;
        else {                    sN = (-d +  b);
                                  sD = a;
        }
    }
    // Finally do the division to get sc and tc:
    sc = std::abs(sN) < SMALL_NUM ? zero : sN / sD;
    tc = std::abs(tN) < SMALL_NUM ? zero : tN / tD;
    
    // Get the difference of the two closest points:
    const auto dP = w + (u * sc) - (v * tc);  // =  S1(sc) - S2(tc)
    
    return dP.dotpr();   // the closest distance
}


template<std::floating_point real> constexpr
auto Geometric<real>::
hexagonal_lattice(
    const A2r orig,
    const real step,
    const szt numLayers
) noexcept -> std::vector<A2r>
{
    std::vector<A2r> v;

    auto add_point = [&](const szt j, const int sign)
    {
        for (szt i=1; i<=2*numLayers+1-j; i++) {
            const auto a1 = half*j + i - numLayers - one;
            const auto a2 = sign * static_cast<real>(j) * std::sin(pi/three);
            v.emplace_back(orig + A2r {a1, a2} * step);
        }
    };

    add_point(0, 0);
    for (szt j=1; j<=numLayers; j++) {
        add_point(j, 1);
        add_point(j, -1);
    }

    return v;
}


template<std::floating_point real> 
constexpr
auto Geometric<real>::
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


template<std::floating_point real> 
constexpr
real Geometric<real>::
grad2rad(
    const real grad
) noexcept
{
    return grad * pi / RAD2GRAD;
}


template<std::floating_point real> 
constexpr
real Geometric<real>::
rad2grad(
    const real rad
) noexcept
{
    return rad * RAD2GRAD / pi;
}


template<std::floating_point real> constexpr
real Geometric<real>::
cos_two_segments(
    const A3r& p1,
    const A3r& p2,
    const A3r& p3,
    const A3r& p4
) noexcept
{
    const auto d1 = p2 - p1;
    const auto d2 = p4 - p3;

    return d1.dotpr(d2) / (d1.norm() * d2.norm());
}


template<std::floating_point real> constexpr
real Geometric<real>::
dotpr_two_segments(
    const A3r& p1,
    const A3r& p2,
    const A3r& p3,
    const A3r& p4
) noexcept
{
    const auto d1 = p2 - p1;
    const auto d2 = p4 - p3;

    return d1.dotpr(d2);
}


// Determines if a point is inside a 2D triangle.
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template<std::floating_point real>
bool Geometric<real>::
point_in_triangle(
    const A2r& p,
    const A2r& v1,
    const A2r& v2,
    const A2r& v3
) noexcept
{
    auto sign = [](const A2r& p1,
                   const A2r& p2,
                   const A2r& p3) noexcept
    {
        const auto d {(p1[0] - p3[0]) * (p2[1] - p3[1]) -
                      (p2[0] - p3[0]) * (p1[1] - p3[1])};
        return d;
    };

    const real d1 {sign(p, v1, v2)};
    const real d2 {sign(p, v2, v3)};
    const real d3 {sign(p, v3, v1)};

    return !((d1 < zero || d2 < zero || d3 < zero) &&
             (d1 > zero || d2 > zero || d3 > zero));
}


// Determines if a point is inside a 2D triangle.
// https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
template<std::floating_point real> 
bool Geometric<real>::
point_in_triangle(
    const real* p,
    const real* v1,
    const real* v2,
    const real* v3
) noexcept
{
    auto sign = [](
        const real* p1,
        const real* p2,
        const real* p3 
    ) noexcept
    {
        return {(*p1 - *p3) * (*(p2+1) - *(p3+1)) -
                (*p2 - *p3) * (*(p1+1) - *(p3+1))};
    };

    const real d1 {sign(p, v1, v2)};
    const real d2 {sign(p, v2, v3)};
    const real d3 {sign(p, v3, v1)};

    return !((d1 < zero || d2 < zero || d3 < zero) &&
             (d1 > zero || d2 > zero || d3 > zero));
}

}  // namespace utils::common

#endif     // UTILS_COMMON_GEOMETRIC_FUNCTIONS
