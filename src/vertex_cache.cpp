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
// mapnik
#include <mapnik/global.hpp>
#include <mapnik/vertex_cache.hpp>
#include <mapnik/offset_converter.hpp>
#include <mapnik/make_unique.hpp>

namespace mapnik
{

vertex_cache::vertex_cache(vertex_cache && rhs)
    : current_position_(std::move(rhs.current_position_)),
      subpaths_(std::move(rhs.subpaths_)),
      position_in_segment_(std::move(rhs.position_in_segment_)),
      angle_(std::move(rhs.angle_)),
      angle_valid_(std::move(rhs.angle_valid_)),
      offseted_lines_(std::move(rhs.offseted_lines_)),
      position_(std::move(rhs.position_))
{
    // The C++11 standard doesn't guarantee iterators are valid when container is moved.
    // We can create them from indexes but we don't need to. Just let them uninitialized.
    initialized_ = false;
}

double vertex_cache::current_segment_angle()
{
    return std::atan2(current_segment_->pos.y - previous_segment_->pos.y,
                      current_segment_->pos.x - previous_segment_->pos.x);
}

double vertex_cache::angle(double width)
{
    double tmp = width + position_in_segment_;
    if ((tmp <= current_segment_->length) && (tmp >= 0))
    {
        //Only calculate angle on request as it is expensive
        if (!angle_valid_)
        {
            angle_ = current_segment_angle();
        }
    }
    else
    {
        scoped_state s(*this);
        if (move(width))
        {
            pixel_position const& old_pos = s.get_state().position();
            return std::atan2(current_position_.y - old_pos.y,
                              current_position_.x - old_pos.x);
        }
        else
        {
            s.restore();
            angle_ = current_segment_angle();
        }
    }
    return width >= 0 ? angle_ : angle_ + M_PI;
}

bool vertex_cache::next_subpath()
{
    if (!initialized_)
    {
        current_subpath_ = subpaths_.begin();
        initialized_ = true;
    }
    else
    {
        current_subpath_++;
    }
    if (current_subpath_ == subpaths_.end()) return false;
    rewind_subpath(); //Initialize position values
    return true;
}

void vertex_cache::rewind_subpath()
{
    current_segment_ = current_subpath_->vector.begin();
    previous_segment_ = current_subpath_->vector.begin();
    //All subpaths contain at least one segment (i.e. the starting point)
    current_position_ = current_segment_->pos;
    position_in_segment_ = 0;
    angle_valid_ = false;
    position_ = 0;
}

void vertex_cache::reset()
{
    initialized_ = false;
}

bool vertex_cache::next_segment()
{
    previous_segment_ = current_segment_; //Next segments starts at the end of the current one
    if (current_segment_ == current_subpath_->vector.end()) return false;
    current_segment_++;
    angle_valid_ = false;
    if (current_segment_ == current_subpath_->vector.end()) return false;
    return true;
}

bool vertex_cache::previous_segment()
{
    if (current_segment_ == current_subpath_->vector.begin()) return false;
    current_segment_--;
    angle_valid_ = false;
    if (current_segment_ == current_subpath_->vector.begin())
    {
        //First segment is special
        previous_segment_ = current_segment_;
        return true;
    }
    previous_segment_--;
    return true;
}

vertex_cache & vertex_cache::get_offseted(double offset)
{
    if (std::fabs(offset) < 0.01)
    {
        return *this;
    }
    offseted_lines_map::iterator pos = offseted_lines_.find(offset);
    if (pos == offseted_lines_.end())
    {
        offset_converter<vertex_cache> converter(*this, false);
        converter.set_offset(offset);
        pos = offseted_lines_.emplace(offset, std::make_unique<vertex_cache>(converter)).first;
    }
    vertex_cache_ptr & offseted_line = pos->second;

    offseted_line->reset();
    offseted_line->next_subpath(); //TODO: Multiple subpath support

    // find the point on the offset line closest to the current position,
    // which we'll use to make the offset line aligned to this one.
    int iter_offset = previous_segment_ - current_subpath_->vector.begin();
    int vec_length = current_subpath_->vector.size();
    offseted_line->position_closest_to(current_position_, current_subpath_->length, position_, iter_offset, vec_length, offset*offset);
    return *offseted_line;
}

inline double dist_sq(pixel_position const &d)
{
    return d.x*d.x + d.y*d.y;
}

void vertex_cache::position_closest_to(pixel_position const &target_pos, 
                                       double unoffset_path_length,
                                       double unoffset_position,
                                       int iter_offset,
                                       int vec_length,
                                       double offset_sq)
{
    // So to speed this up we are making an assumption that the offset line's difference in length
    // compared to the length of the unoffset line gives us a meaningful range to search within.
    //
    // First we find the difference between the two lenghts (this could be negative)
    double length_difference = current_subpath_->length - unoffset_path_length;
    
    // In some cases like a zig zagging line, the difference in total length could be near zero
    // but this will not help find the correct point quickly as it will be shifted in correctly.
    // Therefore we will shift the search range to be a percentage of the total length by
    // overwriting the length difference.
    /*if (length_difference < unoffset_path_length * 0.05)
    {
        length_difference = unoffset_path_length * 0.05; // At most will search 20% of the total length
    }*/

    // The starting position will be based on the expected position and the length difference
    double start_pos = unoffset_position;
    if (length_difference < 0.0)
    {
        start_pos += 2.5 * length_difference;
    }
    else 
    {
        start_pos -= 1.5 * length_difference;
    }
    
    if (start_pos < 0.0)
    {
        start_pos = 0.0;
    }

    int vec_diff = current_subpath_->vector.size() - vec_length;

    if (vec_diff < 0)
    {
        iter_offset += vec_diff;
        if (iter_offset < 0)
        {
            iter_offset = 0;
        }
    }

    // Where to stop searching.
    double end_pos = start_pos + 2.0 * std::abs(length_difference);
    double end_pos2 = start_pos + 3.0 * std::abs(length_difference);
    
    // First subpath is intial point and is special so lets skip forward!
    current_segment_ += iter_offset;
    if (current_segment_ == current_subpath_->vector.begin())
    {
        previous_segment_ = current_segment_;
        current_segment_++;
    }
    else
    {
        previous_segment_ = current_segment_ - 1;
    }
    position_ = previous_segment_->total_length;

    if (start_pos < position_)
    {
        move(start_pos - position_);
    }

    // Save the initial state.
    current_position_ = previous_segment_->pos;
    state s = save_state();
    position_ -= position_in_segment_;
    position_in_segment_ = 0.0; 

    // Finding the closet distance from the segment to the 
    // Find the fartherest position we will search
    double min_dist_sq = std::numeric_limits<double>::max();
    bool last = false;

    // find closest approach of each individual segment to the
    // target position. would be good if there were some kind
    // of prior, or fast test to avoid calculating on each
    // segment, but i can't think of one.
    offset_sq *= 1.2;
    int count = 0;
    do
    {
        count++;
        if (position_ > end_pos && (position_ > end_pos2 || min_dist_sq <= offset_sq))
        {
            last = true;
        }
        
        double t_minus_p_x = target_pos.x - previous_segment_->pos.x;
        double t_minus_p_y = target_pos.y - previous_segment_->pos.y;
        double end_dist_sq = t_minus_p_x*t_minus_p_x + t_minus_p_y*t_minus_p_y;
        if (end_dist_sq < min_dist_sq)
        {
            min_dist_sq = end_dist_sq;
            current_position_ = previous_segment_->pos;
            s = save_state();
        }

        double dx = current_segment_->pos.x - previous_segment_->pos.x;
        double dy = current_segment_->pos.y - previous_segment_->pos.y;
        if ((dx != 0.0) || (dy != 0.0))
        {
            double t = (t_minus_p_x * dx + t_minus_p_y * dy) / (dx*dx + dy*dy);

            if ((t >= 0.0) && (t <= 1.0))
            {
                double pt_x = target_pos.x - ((dx * t) + previous_segment_->pos.x);
                double pt_y = target_pos.y - ((dy * t) + previous_segment_->pos.y);
                double pt_dist_sq = pt_x*pt_x + pt_y*pt_y;

                if (pt_dist_sq < min_dist_sq)
                {
                    min_dist_sq = pt_dist_sq;
                    position_in_segment_ = current_segment_->length * t;
                    position_ += position_in_segment_;
                    double factor = position_in_segment_ / current_segment_->length;
                    current_position_ = previous_segment_->pos + (current_segment_->pos - previous_segment_->pos) * factor;
                    s = save_state();
                    position_ -= position_in_segment_;
                    position_in_segment_ = 0.0; 
                }
            }
        }

        position_ += current_segment_->length;
        previous_segment_ = current_segment_; //Next segments starts at the end of the current one
        current_segment_++;
    }
    while (current_segment_ != current_subpath_->vector.end() && !last);
    restore_state(s);
}

bool vertex_cache::forward(double length)
{
    if (length < 0)
    {
        MAPNIK_LOG_ERROR(vertex_cache) << "vertex_cache::forward() called with negative argument!\n";
        return false;
    }
    return move(length);
}

bool vertex_cache::backward(double length)
{
    if (length < 0)
    {
        MAPNIK_LOG_ERROR(vertex_cache) << "vertex_cache::backward() called with negative argument!\n";
        return false;
    }
    return move(-length);
}

bool vertex_cache::move(double length)
{
    if (current_segment_ == current_subpath_->vector.end()) return false;

    position_ += length;
    length += position_in_segment_;
    while (length >= current_segment_->length)
    {
        length -= current_segment_->length;
        if (!next_segment()) return false; //Skip all complete segments
    }
    while (length < 0)
    {
        if (!previous_segment()) return false;
        length += current_segment_->length;
    }
    double factor = length / current_segment_->length;
    position_in_segment_ = length;
    current_position_ = previous_segment_->pos + (current_segment_->pos - previous_segment_->pos) * factor;
    return true;
}

bool vertex_cache::move_to_distance(double distance)
{
    if (current_segment_ == current_subpath_->vector.end()) return false;

    double position_in_segment = position_in_segment_ + distance;
    if (position_in_segment < .0 || position_in_segment >= current_segment_->length)
    {
        // If there isn't enough distance left on this segment
        // then we need to search until we find the line segment that ends further than distance away
        double abs_distance = std::abs(distance);
        double new_abs_distance = .0;
        pixel_position inner_pos; // Inside circle.
        pixel_position outer_pos; // Outside circle.

        position_ -= position_in_segment_;

        if (distance > .0)
        {
            do
            {
                position_ += current_segment_->length;
                if (!next_segment()) return false;
                new_abs_distance = (current_position_ - current_segment_->pos).length();
            }
            while (new_abs_distance < abs_distance);

            inner_pos = previous_segment_->pos;
            outer_pos = current_segment_->pos;
        }
        else
        {
            do
            {
                if (!previous_segment()) return false;
                position_ -= current_segment_->length;
                new_abs_distance = (current_position_ - previous_segment_->pos).length();
            }
            while (new_abs_distance < abs_distance);

            inner_pos = current_segment_->pos;
            outer_pos = previous_segment_->pos;
        }

        find_line_circle_intersection(current_position_.x, current_position_.y, abs_distance,
            inner_pos.x, inner_pos.y, outer_pos.x, outer_pos.y,
            current_position_.x, current_position_.y);

        position_in_segment_ = (current_position_ - previous_segment_->pos).length();
        position_ += position_in_segment_;
    }
    else
    {
        position_ += distance;
        distance += position_in_segment_;
        double factor = distance / current_segment_->length;
        position_in_segment_ = distance;
        current_position_ = previous_segment_->pos + (current_segment_->pos - previous_segment_->pos) * factor;
    }
    return true;
}

void vertex_cache::rewind(unsigned)
{
    vertex_subpath_ = subpaths_.begin();
    vertex_segment_ = vertex_subpath_->vector.begin();
}

unsigned vertex_cache::vertex(double *x, double *y)
{
    if (vertex_segment_ == vertex_subpath_->vector.end())
    {
        vertex_subpath_++;
        if (vertex_subpath_ == subpaths_.end()) return agg::path_cmd_stop;
        vertex_segment_ = vertex_subpath_->vector.begin();
    }
    *x = vertex_segment_->pos.x;
    *y = vertex_segment_->pos.y;
    unsigned cmd = (vertex_segment_ == vertex_subpath_->vector.begin()) ? agg::path_cmd_move_to : agg::path_cmd_line_to;
    vertex_segment_++;
    return cmd;
}


vertex_cache::state vertex_cache::save_state() const
{
    state s;
    s.current_segment = current_segment_;
    s.previous_segment = previous_segment_;
    s.position_in_segment = position_in_segment_;
    s.current_position = current_position_;
    s.position_ = position_;
    return s;
}

void vertex_cache::restore_state(state const& s)
{
    current_segment_ = s.current_segment;
    previous_segment_ = s.previous_segment;
    position_in_segment_ = s.position_in_segment;
    current_position_ = s.current_position;
    position_ = s.position_;
    angle_valid_ = false;
}

void vertex_cache::find_line_circle_intersection(
    double cx, double cy, double radius,
    double x1, double y1, double x2, double y2,
    double & ix, double & iy) const
{
    double dx = x2 - x1;
    double dy = y2 - y1;

    double A = dx * dx + dy * dy;
    double B = 2 * (dx * (x1 - cx) + dy * (y1 - cy));
    double C = (x1 - cx) * (x1 - cx) + (y1 - cy) * (y1 - cy) - radius * radius;

    double det = B * B - 4 * A * C;
    if (A <= 1.0e-7 || det < 0)
    {
        // Should never happen.
        // No real solutions.
        return;
    }
    else if (det == 0)
    {
        // Could potentially happen....
        // One solution.
        double t = -B / (2 * A);
        ix = x1 + t * dx;
        iy = y1 + t * dy;
        return;
    }
    else
    {
        // Two solutions.

        // Always use the 1st one
        // We only really have one solution here, as we know the line segment will start in the circle and end outside
        double t = (-B + std::sqrt(det)) / (2 * A);
        ix = x1 + t * dx;
        iy = y1 + t * dy;

        //t = (-B - std::sqrt(det)) / (2 * A);
        //ix = x1 + t * dx;
        //iy = y1 + t * dy;

        return;
    }
}

} //ns mapnik
