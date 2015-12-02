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

#ifndef MAPNIK_OFFSET_CONVERTER_HPP
#define MAPNIK_OFFSET_CONVERTER_HPP

#ifdef MAPNIK_LOG
#include <mapnik/debug.hpp>
#endif
#include <mapnik/global.hpp>
#include <mapnik/config.hpp>
#include <mapnik/vertex.hpp>
#include <mapnik/offset_utils.hpp>
#include <mapnik/geometry_circular_list.hpp>

// boost
#pragma GCC diagnostic push
#include <mapnik/warning_ignore.hpp>
#include <boost/version.hpp>
#pragma GCC diagnostic pop

// stl
#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace mapnik
{

template <typename Geometry>
struct offset_converter
{
    enum status
    {
        initial,
        processed,
        processed_with_error
    };

private:
    Geometry & geom_;
    mapnik::geometry::clist<double> new_geom_;
    double offset_;
    double miter_limit_;
    double arc_tolerance_;
    mapnik::geometry::corner_type corner_;
    status status_;

public:
    offset_converter(Geometry & geom)
        : geom_(geom)
        , offset_(0.0)
        , miter_limit_(mapnik::geometry::DEFAULT_MITER_LIMIT)   
        , arc_tolerance_(mapnik::geometry::DEFAULT_ARC_TOLERANCE)
        , corner_(mapnik::geometry::CORNER_ROUND_TYPE)
        , status_(initial)
    {}

    unsigned type() const
    {
        return static_cast<unsigned>(geom_.type());
    }

    double get_offset() const
    {
        return offset_;
    }

    double get_arc_tolerance() const
    {
        return arc_tolerance_;
    }

    bool get_miter_limit() const
    {
        return miter_limit_;
    }

    mapnik::geometry::corner_type get_corner_type() const
    {
        return corner_;
    }

    void set_offset(double value)
    {
        if (offset_ != value)
        {
            offset_ = value;
            if (status_ != initial)
            {
                reset();
            }
        }
    }

    void set_arc_tolerance(double value)
    {
        if (arc_tolerance_ != value)
        {
            arc_tolerance_ = value;
            if (status_ != initial)
            {
                reset();
            }
        }
    }

    void set_miter_limit(double value)
    {
        if (miter_limit_ != value)
        {
            miter_limit_ = value;
            if (status_ != initial)
            {
                reset();
            }
        }
    }

    void set_corner_type(mapnik::geometry::corner_type value)
    {
        if (corner_ != value)
        {
            corner_ = value;
            if (status_ != initial)
            {
                reset();
            }
        }
    }

    unsigned vertex(double * x, double * y)
    {
        if (offset_ == 0.0)
        {
            return geom_.vertex(x, y);
        }

        if (status_ == initial)
        {
            init_offset_geom();
        }
        // temp for now
        return geom_.vertex(x,y);
    }

    void reset()
    {
        geom_.rewind(0);
        new_geom_.clear();
        status_ = initial;
    }

    void rewind(unsigned)
    {
        // do something here?
    }

private:

    void init_offset_geom()
    {
        // Check if the data has already been intialized so we can skip processing.
        if (status_ == processed || status_ == processed_with_error) 
        {
            return;
        }

        bool success = new_geom_.add_streaming_offset_geometry(
                                 geom_,
                                 offset_,
                                 corner_,
                                 miter_limit_,
                                 arc_tolerance_
                                );
        if (success)
        {
            status_ = processed;
        }
        else
        {
            status_ = processed_with_error;
            return;
        }
    }

};

}

#endif // MAPNIK_OFFSET_CONVERTER_HPP
