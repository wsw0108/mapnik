/*****************************************************************************
 *
 * This file is part of Mapnik (c++ mapping toolkit)
 *
 * Copyright (C) 2015 Artem Pavlenko
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *****************************************************************************/

#ifndef MAPNIK_OFFSET_UTILS_HPP
#define MAPNIK_OFFSET_UTILS_HPP

// Many of the concepts around offseting and the offsetting code have come
// from:

// The Angus Clipping Library
// http://sourceforge.net/projects/polyclipping/files/
// http://www.angusj.com/delphi/clipper.php

// POLYGON OFFSETTING BY COMPUTING WINDING NUMBERS
// http://www.me.berkeley.edu/~mcmains/pubs/DAC05OffsetPolygon.pdf

// std
#include <cmath>

// mapnik
#include <mapnik/vertex.hpp>
#include <mapnik/geometry.hpp>


namespace mapnik
{

namespace geometry
{

enum corner_type : std::uint8_t
{
    CORNER_SQUARE_TYPE = 0,
    CORNER_ROUND_TYPE,
    CORNER_MITER_TYPE,
    CORNER_TYPE_MAX
};

static double DEFAULT_MITER_LIMIT = 2.0;
// See http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Classes/ClipperOffset/Properties/ArcTolerance.htm
// and offset_triginometry2.svg from angus clipper documentation
static double DEFAULT_ARC_TOLERANCE = 0.25;

class offsetter
{
private:
    double miter_limit_;
    // The following is taken from offset_triginometry3.svg in angus clipper documentation
    /*****************************************************
    *
    *                  b (offset point)
    *                  @
    *                 /|\
    *                / | \      ß = angle formed by points (lam)
    *               /  |  \
    *              /   |   @ l
    *             /    |    \
    *            /     |   _ @ c <--- A RIGHT ANGLE
    *           /      | _/   \       pardon the ASCI
    *          /     a @/      \      art.
    *         /       / \       \
    *        /       / Ø \       \    length(ac) == offset
    *       /       /     \       \
    *      /       /       @ m     \
    *******************************************************/
    // When 'mitering' offset polygons, the maximum distance point 'b' can be from 
    // point 'a' is set by 'limit' where limit is a multiple of 'offset'. Therefore, 
    // for any given angle we need to know if length(ab) > limit * offset.
    //
    // Find the largest angle ß (or smallest Ø since Ø = pi - ß) for a given 
    // limit, expressing ß as sin(ß) or cos(ß) since these can easily be derived 
    // from cross or dot products respectively. 
    //
    // angle(abc) = Ø/2 
    // length(ab) = limit * offset
    // length(ac) = offset
    // sin(Ø/2) = offset / (limit * offset) = 1 / limit
    // Given that sin(Ø/2) = sqrt((1-cos(Ø))/2) **
    // 1 / limit = sqrt((1-cos(Ø))/2)
    // limit = sqrt(2 / (1-cos(Ø)))
    // 1-cos(Ø) = 2 / sqr(limit) 
    // Since Ø = pi - ß ...
    // 1-cos(pi - ß) = 2 / sqr(limit)
    // and given cos(pi-ß) = -cos(ß) ** ... 
    // 1+cos(ß) = 2 / sqr(limit) 
    // cos(ß) = 2 / sqr(limit) - 1
    
    // Example: if miter limit = 2 (ie 2 times offset) then cos(ß) = (2 / 4) -1 = -0.5 
    // and so ß = 120 degrees. Therefore, when ß > 120 deg. (or Ø < 60 deg.), the 
    // distance point 'b' would be from point 'a' would exceed the limit.
    
    double m_limit_sqr_over_two_; // This is 2 / sqr(limit)
    double arc_tolerance_;
    double steps_sin_;
    double steps_cos_;
    double steps_per_radian_;
    double offset_;
    corner_type corner_;

public:
    offsetter(double offset,
             corner_type corner = CORNER_SQUARE_TYPE,
             double miter_limit = DEFAULT_MITER_LIMIT,   
             double arc_tolerance = DEFAULT_ARC_TOLERANCE)
        : miter_limit_(miter_limit),
          m_limit_sqr_over_two_(0.5),
          arc_tolerance_(arc_tolerance),
          steps_sin_(0.0),
          steps_cos_(0.0),
          steps_per_radian_(0.0),
          offset_(offset),
          corner_(corner)
    {
        update_miter_limit();
        update_arc_tolerance();
    }

    template <typename Vtx>
    static point<double> unit_normal(Vtx const& v0, Vtx const& v1)
    {
        // Determine dy and dx
        double dx = v0.x - v1.x; 
        double dy = v0.y - v1.y;
        double ds = std::sqrt(dx*dx + dy*dy); // magnitude of segment vector
        dx = dx / ds;
        dy = dy / ds;
        // Unit normal vector to the LEFT.
        return point<double>(-dy, dx);
    }

    template <typename Vtx>
    point<typename Vtx::value_type> offset_dx_dy(Vtx const& v, double x, double y)
    {
        using value_type = typename Vtx::value_type;
        return point<value_type>(v.x + static_cast<value_type>(std::round(x)), 
                                 v.y + static_cast<value_type>(std::round(y)));
    }

    point<double> offset_dx_dy(point<double> const& v, double x, double y)
    {
        return point<double>(v.x + x, v.y + y);
    }

    point<double> offset_dx_dy(vertex2d const& v, double x, double y)
    {
        return point<double>(v.x + x, v.y + y);
    }
    
    template <typename Vtx>
    point<typename Vtx::value_type> offset_point(Vtx const& v0,
                                        point<double> const& seg0_unit_normal)
    {
        // Solve for the point at p0 that is offset from
        // the point at v0. The point is offset by the distance
        // "offset" on the unit vector normal to the segment seg0.
        /*******************************************************
        *
        *               o 
        *              /
        *             /
        *    p0      / seg0
        *     @     /
        *          o
        *          v0
        * 
        * Note: showing a positive offset here.
        *******************************************************/
        
        return offset_dx_dy(v0, 
                            offset_ * seg0_unit_normal.x, 
                            offset_ * seg0_unit_normal.y);
    }
    
    template <typename Vtx, typename Container>
    void offset_corner(Vtx const& v1,
                       point<double> const& seg0_unit_normal,
                       point<double> const& seg1_unit_normal,
                       Container & container)
    {
        /***********************************************
        * We are creating the offset section around the 
        * corner of a joint of two segments as shown below:
        *            
        *                 v1    seg1   v2
        *                o--------------o
        *               / 
        *              /               
        *             / seg0               
        *            /                 
        *           o                 
        *           v0                 
        *
        * The side of the offset segments required will
        * will be based on the offset being negative
        * or positive.
        ************************************************/
        
        // We can find the angle between the two unit normals
        // relative to the v1 in order to solve for the angle
        // and offset set of points. We can call the angle
        // between the two unit vectors to be "theta".

        // The determinate of the two unit normal vectors is equal 
        // to "magnitude" of sin(Ø).
        double m_sin_theta = seg0_unit_normal.x * seg1_unit_normal.y - seg0_unit_normal.y * seg1_unit_normal.x;
        double m_cos_theta = seg0_unit_normal.x * seg1_unit_normal.x + seg0_unit_normal.y * seg1_unit_normal.y;
        
        // If absolute value of m_sin_theta multiplied by offset is less then 1, we know that 
        // the hypotenuse calculated using the two angles will be less then the value of our offset
        if (std::fabs(m_sin_theta * offset_) < 1.0)
        {
            // Now we need to check that the angle is very near to either 2pi or 0
            // and not on the other side of the unit circle. Therefore, we need
            // to check if cos would put it in quadrants II or III. 
            if (m_cos_theta > 0.0)
            {
                // We know this is in the I or IV quadrant, therefore this is a very slight
                // angle. Therefore we should simply just offset from unit normal of seg1.
                // (seg0 would work as well more then likely)
                // This is basically a straight line
                container.push_back(offset_point(v1, seg1_unit_normal));
                return;
            }
        }
        // It seems confusing to cap sin at this
        // point but the purpose is that we are
        // now calculating ß rather then Ø. This
        // will be used in the case we have 
        // a square or round corner.
        else if (m_sin_theta > 1.0) 
        {
            m_sin_theta = 1.0;
        }
        else if (m_sin_theta < -1.0)
        {
            m_sin_theta = -1.0;
        }

        // Next we need to check if the angle would be concave
        // The previous capping of sine won't change the result of this
        if (m_sin_theta * offset_ < 0.0)
        {
            // If the angle is concave we add three offset points
            // it is expected that they will cause intersections
            // but a winding order calculation later will fix these
            // sort of issues.
            offset_corner_concave(v1, seg0_unit_normal, seg1_unit_normal, container);
            return;
        }
        
        // The angle is going to be convex so we can decide the type of corner to apply
        double cos_beta_plus_1 = 1 + m_cos_theta;
        switch (corner_)
        {
        case CORNER_SQUARE_TYPE:
            offset_corner_square(m_sin_theta,
                                 m_cos_theta,
                                 v1,
                                 seg0_unit_normal,
                                 seg1_unit_normal,
                                 container);
            break;
        case CORNER_MITER_TYPE:
            // See definition of m_limit_sqr_over_two_ above
            if (cos_beta_plus_1 >= m_limit_sqr_over_two_)
            {
                offset_corner_miter(cos_beta_plus_1,
                                    v1,
                                    seg0_unit_normal,
                                    seg1_unit_normal,
                                    container);
            }
            else
            {
                offset_corner_square(m_sin_theta,
                                     m_cos_theta,
                                     v1,
                                     seg0_unit_normal,
                                     seg1_unit_normal,
                                     container);
            }
            break;
        case CORNER_TYPE_MAX:
        case CORNER_ROUND_TYPE:
        default:
            offset_corner_round(m_sin_theta,
                                m_cos_theta,
                                v1,
                                seg0_unit_normal,
                                seg1_unit_normal,
                                container);
            break;
        }
    }
    
    template <typename Vtx, typename Container>
    void offset_corner_concave(Vtx const& v1,
                               point<double> const& seg0_unit_normal,
                               point<double> const& seg1_unit_normal,
                               Container & container)
    {
        container.push_back(offset_point(v1, seg0_unit_normal));
        container.push_back(point<typename Vtx::value_type>(v1.x, v1.y));
        container.push_back(offset_point(v1, seg1_unit_normal));
    }
    
    template <typename Vtx, typename Container>
    void offset_corner_square(double m_sin_theta,
                              double m_cos_theta,
                              Vtx const& v1,
                              point<double> const& seg0_unit_normal,
                              point<double> const& seg1_unit_normal,
                              Container & container)
    {
        // We want to calculate tan(ß/4) as this gives us the distance
        // to the point where the squaring will start.
        // See angus clipper documentation offset_triginometry.svg
        
        double beta = std::atan2(m_sin_theta, m_cos_theta);
        double tan_beta_4 = std::tan(beta / 4);
        
        container.push_back(offset_dx_dy(
                                v1,
                                offset_ * (seg0_unit_normal.x - seg0_unit_normal.y * tan_beta_4),
                                offset_ * (seg0_unit_normal.y + seg0_unit_normal.x * tan_beta_4)
                            ));
        container.push_back(offset_dx_dy(
                                v1,
                                offset_ * (seg1_unit_normal.x + seg1_unit_normal.y * tan_beta_4),
                                offset_ * (seg1_unit_normal.y - seg1_unit_normal.x * tan_beta_4)
                            ));
    }

    template <typename Vtx, typename Container>
    void offset_corner_miter(double cos_beta_plus_one,
                             Vtx const& v1,
                             point<double> const& seg0_unit_normal,
                             point<double> const& seg1_unit_normal,
                             Container & container)
    {
        double q = offset_ / cos_beta_plus_one;
        container.push_back(offset_dx_dy(
                                v1,
                                q * (seg0_unit_normal.x + seg1_unit_normal.x),
                                q * (seg0_unit_normal.y + seg1_unit_normal.y)
                            ));
    }

    template <typename Vtx, typename Container>
    void offset_corner_round(double m_sin_theta,
                             double m_cos_theta,
                             Vtx const& v1,
                             point<double> const& seg0_unit_normal,
                             point<double> const& seg1_unit_normal,
                             Container & container)
    {
        double beta = std::atan2(m_sin_theta, m_cos_theta);
        std::size_t steps = static_cast<std::size_t>(std::round(steps_per_radian_ * std::fabs(beta)));
        
        // Always have at least one step
        if (steps < 1)
        {
            steps = 1;
        }
        
        // Vector unit vector components to the position around the circle
        // starting at the unit normal vector of seg0.
        double vec_x = seg0_unit_normal.x;
        double vec_y = seg0_unit_normal.y;

        // Add the first position.
        container.push_back(offset_point(v1, seg0_unit_normal));

        for (std::size_t i = 1; i < steps; ++i)
        {
            // Calculate next vector position
            double prev_vec_x = vec_x;
            vec_x = prev_vec_x * steps_cos_ - vec_y * steps_sin_;
            vec_y = prev_vec_x * steps_sin_ + vec_y * steps_cos_;
            container.push_back(offset_dx_dy(
                                    v1,
                                    offset_ * vec_x,
                                    offset_ * vec_y
                                ));
        }

        // Now add the last point which is an offset from seg1.
        container.push_back(offset_point(v1, seg1_unit_normal));
    }

    template <typename Vtx, typename Container>
    void offset_point_circle(Vtx const& v1, Container & container)
    {
        std::size_t steps = static_cast<std::size_t>(std::round(steps_per_radian_ * 2.0 * M_PI));
        
        // Vector unit vector components to the position around the circle
        // starting at the unit normal vector of seg0.
        double vec_x = 1;
        double vec_y = 0;

        // Add the first position.
        container.push_back(offset_dx_dy(
                                v1, 
                                offset_ * vec_x,
                                offset_ * vec_y
                            ));
        
        // Finish out the circle
        for (std::size_t i = 1; i < steps; ++i)
        {
            // Calculate next vector position
            double prev_vec_x = vec_x;
            vec_x = prev_vec_x * steps_cos_ - vec_y * steps_sin_;
            vec_y = prev_vec_x * steps_sin_ + vec_y * steps_cos_;
            container.push_back(offset_dx_dy(
                                    v1,
                                    offset_ * vec_x,
                                    offset_ * vec_y
                                ));
        }

    }

private:

    void update_miter_limit()
    {
        if (miter_limit_ > 2.0)
        {
            m_limit_sqr_over_two_ = 2.0 / (miter_limit_ * miter_limit_);
        }
        else
        {
            m_limit_sqr_over_two_ = 0.5;
        }
    }

    void update_arc_tolerance()
    {
        double arc_tol;
        if (arc_tolerance_ <= 0.0) 
        {
            arc_tol = DEFAULT_ARC_TOLERANCE;
        }
        else if (arc_tolerance_ > std::fabs(offset_) * DEFAULT_ARC_TOLERANCE) 
        {
            arc_tol = std::fabs(offset_) * DEFAULT_ARC_TOLERANCE;
        }
        else 
        {
            arc_tol = arc_tolerance_;
        }
        double steps = M_PI / std::acos(1 - arc_tol / std::fabs(offset_));
        if (steps > std::fabs(offset_) * M_PI)
        {
            steps = std::fabs(offset_) * M_PI;  //ie excessive precision check
        }
        steps_sin_ = std::sin((2.0 * M_PI) / steps);
        steps_cos_ = std::cos((2.0 * M_PI) / steps);
        steps_per_radian_ = steps / (2.0 * M_PI);
        if (offset_ < 0.0)
        {
            steps_sin_ = -steps_sin_;
        }
    }
};

} // end ns geometry

} // end ns mapnik

#endif // MAPNIK_OFFSET_UTILS_HPP
