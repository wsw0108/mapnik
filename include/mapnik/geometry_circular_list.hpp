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

#ifndef MAPNIK_GEOMETRY_CIRCULAR_LIST_HPP
#define MAPNIK_GEOMETRY_CIRCULAR_LIST_HPP

// mapnik
#include <mapnik/util/variant.hpp>
#include <mapnik/geometry.hpp>
// boost
#include <boost/intrusive/circular_list_algorithms.hpp>

namespace mapnik
{

namespace geometry
{

enum clist_type : std::uint8_t
{
    CLIST_POINT_TYPE = 0,
    CLIST_LINE_TYPE,
    CLIST_RING_TYPE,
    CLIST_NULL_TYPE,
    CLIST_TYPE_MAX
};

template <typename T>
struct clist_vertex
{
    clist_vertex(T x_,
                 T y_,
                 std::size_t part_index_)
        : p(x_, y_),
          part_index(part_index_),
          next_(nullptr),
          prev_(nullptr) {}
    
    clist_vertex(point<T> p_,
                 std::size_t part_index_)
        : p(p_),
          part_index(part_index_),
          next_(nullptr),
          prev_(nullptr) {}
    
    clist_vertex(point<T> && p_,
                 std::size_t part_index_)
        : p(p_),
          part_index(part_index_),
          next_(nullptr),
          prev_(nullptr) {}

    bool operator <(clist_vertex<T> const& v) const
    {
        return p.x < v.p.x;
    }

    point<T> p;
    std::size_t part_index;
    clist_vertex<T> * next_;
    clist_vertex<T> * prev_;
};

template <typename T>
struct clist_vertex_node_traits
{
    typedef clist_vertex<T> node;
    typedef clist_vertex<T> * node_ptr;
    typedef const clist_vertex<T> * const_node_ptr;

    static node_ptr get_next(const_node_ptr n)         
    {  
        return n->next_;  
    }

    static void set_next(node_ptr n, node_ptr next)    
    {  
        n->next_ = next;  
    }
    
    static node *get_previous(const_node_ptr n)
    {  
        return n->prev_;  
    }

    static void set_previous(node_ptr n, node_ptr prev)
    {  
        n->prev_ = prev;  
    }
};


template <typename T>
struct clist_ring
{
    clist_ring():
        winding_number(0),
        origin(nullptr) {}
    
    clist_ring(clist_vertex<T> * vertex):
        winding_number(0),
        origin(vertex) {}
    
    int winding_number;
    // A link to a vertex on the ring that will be its origin.
    clist_vertex<T> * origin;

    bool empty() const
    {
        return origin == nullptr;
    }
};

template <typename T>
struct clist_line
{
    using clist_algo = boost::intrusive::circular_list_algorithms<clist_vertex_node_traits<T> >;
    
    clist_line():
        origin_(0,0,0) {}
    
    clist_line(clist_vertex<T> * vertex):
        origin_(0,0,0)
    {
        clist_algo::link_before(vertex, &origin_);
    }
    
    bool empty() const
    {
        return origin_ == origin_->next_;
    }

    // A line will be a circular list with the origin as an imaginary
    // vertex. This vertex will not be in the collection's vertices list
    // and is used to link the begining and end points as this is still
    // a circular list
    clist_vertex<T> origin_;

};

template <typename T>
struct clist_point
{
    clist_point():
        origin(nullptr) {}

    clist_point(clist_vertex<T> * vertex):
        origin(vertex) {}
     
    // A link to a vertex on the ring that will be its origin.
    clist_vertex<T> * origin;
       
    bool empty() const
    {
        return origin == nullptr;
    }
};

struct clist_null
{
    bool empty() const
    {
        return true;
    }
};

template <typename T>
using clist_geometry_base = mapnik::util::variant<clist_null,
                                                  clist_point<T>,
                                                  clist_line<T>,
                                                  clist_ring<T> >;
namespace detail
{

struct get_clist_part_empty_visitor
{
    template <typename T>
    bool operator()(T const& data) const
    {
        return data.empty();
    }
};

} // end ns detail

template <typename T>
struct clist_geometry;

template <typename T>
clist_geometry<T> make_clist_geometry_by_type(clist_vertex<T> * vertex, clist_type type)
{
    switch (type)
    {
    default:
    case CLIST_TYPE_MAX:
    case CLIST_NULL_TYPE:
        return clist_null();
    case CLIST_POINT_TYPE:
        return clist_point<T>(vertex);
    case CLIST_LINE_TYPE:
        return clist_line<T>(vertex);
    case CLIST_RING_TYPE:
        return clist_ring<T>(vertex);
    }
}


template <typename T>
struct clist_geometry : clist_geometry_base<T>
{
    using value_type = T;
    
    clist_geometry()
        : clist_geometry_base<value_type>() {} // null geometry

    clist_geometry(clist_vertex<value_type> * vertex, clist_type type)
        : clist_geometry_base<value_type>(std::move(make_clist_geometry_by_type<value_type>(vertex, type))) {}

    template <typename G>
    clist_geometry(G && geom)
        : clist_geometry_base<value_type>(std::forward<G>(geom)) {}
    
    bool empty() const
    {
        return util::apply_visitor(detail::get_clist_part_empty_visitor(),*this);
    }
};

template <typename T>
class clist_stream_container
{
    using value_type = T;
    using clist_algo = boost::intrusive::circular_list_algorithms<clist_vertex_node_traits<value_type> >;
    using clist_vertex_t = clist_vertex<value_type>;
    using clist_vertex_ptr = clist_vertex<value_type> *;
    using point_t = point<value_type>;
public:

    clist_stream_container(std::vector<clist_vertex_t> & vertices, 
                           std::size_t id)
        :  vertices_(vertices),
           id_(id),
           back_(nullptr),
           first_(nullptr) {}

    void push_back(point_t && p)
    {
        if (back_ != nullptr)
        {
            vertices_.emplace_back(p, id_);
            clist_vertex_ptr next = &(vertices_.back());
            clist_algo::link_after(back_, next);
            back_ = next;
        }
        else
        {
            vertices_.emplace_back(p, id_);
            first_ = &(vertices_.back());
            back_ = first_;
            clist_algo::init_header(first_);
        }
    }

    clist_vertex_ptr begin()
    {
        return first_;
    }

    std::size_t get_id()
    {
        return id_;
    }

private:
    std::vector<clist_vertex_t> & vertices_;
    std::size_t id_;
    clist_vertex_ptr back_;
    clist_vertex_ptr first_;
};

template <typename T>
struct clist
{
    using value_type = T;
    using clist_vertex_t = clist_vertex<value_type>;
    using clist_vertex_ptr = clist_vertex<value_type> *;
    using clist_geometry_t = clist_geometry<value_type>;
    using point_t = point<value_type>;
    using container_t = clist_stream_container<value_type>;
    using container_ptr = std::shared_ptr<container_t>;
    using clist_algo = boost::intrusive::circular_list_algorithms<clist_vertex_node_traits<value_type> >;
private:
    std::vector<clist_vertex_t> vertices_;
    std::vector<clist_geometry_t> geometries_;
public:

    clist():
        vertices_(),
        geometries_() {}


    
    // returns TRUE on success and FALSE on a failure
    // If a failure occurs the clist should be considered invalid.
    template <typename Geom>
    bool add_streaming_geometry(Geom & geometry)
    {
        // Position of most current vertex
        vertex<value_type,2> v0(vertex<value_type,2>::no_init);
        // Position of current vertex - 1 in polygon part, if this doesn't exist will be first point in part
        vertex<value_type,2> v1(vertex<value_type,2>::no_init);
        
        /***********************************************
        *
        *             v1    seg0     v0
        *              o-------------o
        *
        ************************************************/

        // Find the first MOVETO prior to processing
        while ((v0.cmd = geometry.vertex(&v0.x, &v0.y)) != SEG_MOVETO)
        {
            if (v0.cmd == SEG_END)
            {
                // If we found segment end here it means that 
                // we found a SEG_END prior to the first moveto
                // so we will just return as this could be empty geometry
                return false;
            }
        }
        
        // Sometimes we might throw out a vertex and "skip" it and not
        // create a segment with it.
        bool skipped_vertex = false;

        // We want to keep track of the first vertex we create often
        // so that a container can link to it.
        std::size_t geometry_id = create_null();
        clist_vertex_ptr first = add_first_vertex(v0.x, v0.y, geometry_id);
        clist_vertex_ptr previous = first;

        for (;;)
        {
            if (!skipped_vertex)
            {
                v1 = v0;
            }
            else
            {
                skipped_vertex = false;
            }
            v0.cmd = geometry.vertex(&v0.x, &v0.y);
            if (v0.cmd == SEG_MOVETO || v0.cmd == SEG_END)
            {
                if (v1.cmd == SEG_LINETO)
                {
                    // mark previous geometry as line
                    geometry_is_line(first, geometry_id);
                }
                else if (v1.cmd == SEG_MOVETO)
                {
                    // mark previous geometry as point
                    geometry_is_point(first, geometry_id);
                    
                }
                if (v0.cmd == SEG_END)
                {
                    break;
                }
                else // v0.cmd == SEG_MOVETO
                {
                    geometry_id = create_null();
                    first = add_first_vertex(v0.x, v0.y, geometry_id);
                    previous = first;
                }
            }
            else if (v0.cmd == SEG_LINETO)
            {
                if (!(v1.cmd == SEG_MOVETO || v1.cmd == SEG_LINETO))
                {
                    // skip vertex because command is invalid
                    skipped_vertex = true;
                }
                else if (v1.x == v0.x && v1.y == v0.y)
                {
                    // line has no length skip it.
                    skipped_vertex = true;
                }
                else 
                {
                    previous = add_vertex(v0.x, v0.y, geometry_id, previous);
                }
            }
            else if (v0.cmd == SEG_CLOSE)
            {
                if (v1.cmd != SEG_LINETO)
                {
                    // Could possibly be two closes in a row which would be
                    // invalid. Could also be a MOVETO followed by a CLOSE.
                    // this would simply be a point so still invalid and
                    // worth skipping.
                    skipped_vertex = true;
                }
                else
                {
                    // If the previous point is already equal
                    // to the first we need to remove it.
                    if (v1.x == first->p.x && v1.y == first->p.y)
                    {
                        previous = remove_last_vertex();
                    }
                    geometry_is_ring(first, geometry_id);
                }
            }
            else
            {
                // Invalid command
                return false;
            }
        }
        return true;
    }

    // returns TRUE on success and FALSE on a failure
    // If a failure occurs the clist should be considered invalid.
    template <typename Geom>
    bool add_streaming_offset_geometry(Geom & geometry,
                                       double offset_distance,
                                       corner_type corner = CORNER_ROUND_TYPE,
                                       double miter_limit = DEFAULT_MITER_LIMIT,   
                                       double arc_tolerance = DEFAULT_ARC_TOLERANCE)
    {
        if (offset_distance == 0.0)
        {
            // if the offset is zero, we do not need to have complex logic to
            // calculate offsets therefore we just use the standard steam into
            // clist logic.
            return add_streaming_geometry(geometry);
        }

        // Now we need to construct a vertex offseter.
        offsetter offset(offset_distance, corner, miter_limit, arc_tolerance); 

        // Position of most current vertex
        vertex<value_type,2> v0(vertex<value_type,2>::no_init);
        // Position of current vertex - 1 in polygon part, if this doesn't exist will be first point in part
        vertex<value_type,2> v1(vertex<value_type,2>::no_init);
        // Position of current vertex - 2 in polygon part, if this doesn't exist will be first point in part
        vertex<value_type,2> v2(vertex<value_type,2>::no_init);
        
        /***********************************************
        *
        *             v1      seg1   v0
        *              o-------------o
        *             /               
        *      seg0  /                  
        *           /                  
        *          o                   
        *         v2                   
        *
        * PLEASE NOTE THAT THE ORDER OF SEGMENTS (v2,v1,v0)
        * IN THE CLIST STREAM CODE IS DIFFERENT THEN
        * THAT IN THE OFFSET CODE (v0,v1,v2)
        ************************************************/

        // Find the first MOVETO prior to processing
        while ((v0.cmd = geometry.vertex(&v0.x, &v0.y)) != SEG_MOVETO)
        {
            if (v0.cmd == SEG_END)
            {
                // If we found segment end here it means that 
                // we found a SEG_END prior to the first moveto
                // so we will just return as this could be empty geometry
                return false;
            }
        }
        
        // Sometimes we might throw out a vertex and "skip" it and not
        // create a segment with it.
        bool skipped_vertex = false;
        
        point<double> seg0_unit_normal(0,0);
        point<double> seg1_unit_normal(0,0);

        container_ptr current_geom;
        
        vertex<value_type,2> first_vertex = v0;
        vertex<value_type,2> second_vertex = v0;

        for (;;)
        {
            if (!skipped_vertex)
            {   
                v2 = v1;
                v1 = v0;
                seg0_unit_normal = seg1_unit_normal;
            }
            else
            {
                skipped_vertex = false;
            }
            v0.cmd = geometry.vertex(&v0.x, &v0.y);
            if (v0.cmd == SEG_MOVETO || v0.cmd == SEG_END)
            {
                // MOVETO and END are very similar in logic
                // the only difference is that END results in a break
                if (v1.cmd == SEG_LINETO)
                {
                    if (v2.cmd == SEG_MOVETO)
                    {
                        // No points have yet been added to this
                        // offset. So we need to start the geometry and add new points.
                        seg0_unit_normal = offsetter::unit_normal(v2,v1);
                        current_geom->push_back(offset.offset_point(v2, seg0_unit_normal));
                        current_geom->push_back(offset.offset_point(v1, seg0_unit_normal));
                    }
                    else // v2.cmd == SEG_LINETO
                    {
                        // We haven't added the segment for the last point either
                        // so we need to add that first.
                        seg0_unit_normal = offsetter::unit_normal(v2, v1);
                        current_geom->push_back(offset.offset_point(v1, seg0_unit_normal));
                        // Because we were concerned this might be a polygon
                        // we didn't calculate the offset of the first vertex
                        // of the line so we need to do it now.
                        seg1_unit_normal = offsetter::unit_normal(first_vertex, second_vertex);
                        current_geom->push_back(offset.offset_point(first_vertex, seg1_unit_normal));
                    }
                    // mark previous geometry as line
                    geometry_is_line(current_geom);
                }
                else if (v1.cmd == SEG_MOVETO)
                {
                    // No need to offset just add the point
                    current_geom->push_back(point_t(v1.x, v1.y));
                    // mark previous geometry as point
                    geometry_is_point(current_geom);
                }
                else if (v1.cmd == SEG_CLOSE)
                {
                    geometry_is_ring(current_geom);
                }

                if (v0.cmd == SEG_END)
                {
                    break;
                }
                else // v0.cmd == SEG_MOVETO
                {
                    current_geom = new_geometry();
                    first_vertex = v0;
                }
            }
            else if (v0.cmd == SEG_LINETO)
            {
                if (!(v1.cmd == SEG_MOVETO || v1.cmd == SEG_LINETO))
                {
                    // skip vertex because command is invalid
                    skipped_vertex = true;
                    continue;
                }
                else if (v1.x == v0.x && v1.y == v0.y)
                {
                    // line has no length skip it.
                    skipped_vertex = true;
                    continue;
                }
                seg1_unit_normal = offsetter::unit_normal(v1,v0);
                if (v1.cmd == SEG_MOVETO)
                {
                    // This only means we have 1 segment.
                    // because we are not sure if this is a ring
                    // or line we will hold off putting an offset point
                    // on the first vertex for now.
                    second_vertex = v0;
                }
                else // v1.cmd == SEG_LINETO
                {
                    offset.offset_corner(v1, seg0_unit_normal, seg1_unit_normal, *current_geom);
                }
            }
            else if (v0.cmd == SEG_CLOSE)
            {
                if (v1.cmd != SEG_LINETO)
                {
                    // Could possibly be two closes in a row which would be
                    // invalid. Could also be a MOVETO followed by a CLOSE.
                    // this would simply be a point so still invalid and
                    // worth skipping.
                    skipped_vertex = true;
                }
                else if (v2.cmd == SEG_MOVETO)
                {
                    // This is a situation with a MOVETO,LINETO,CLOSE
                    // We will just offset the line.
                    seg0_unit_normal = offsetter::unit_normal(v1,v0);
                    seg1_unit_normal = offsetter::unit_normal(v0,v1);
                    offset.offset_corner(v1, seg0_unit_normal, seg1_unit_normal, *current_geom);
                    offset.offset_corner(v0, seg1_unit_normal, seg0_unit_normal, *current_geom);
                }
                else // v2.cmd == SEG_LINETO
                {
                    // If the previous point is already equal
                    // to the first point in the geometry
                    // we need to use slightly different logic
                    if (v1.x == first_vertex.x && v1.y == first_vertex.y)
                    {
                        seg1_unit_normal = offsetter::unit_normal(v1,second_vertex);
                        offset.offset_corner(v1, seg0_unit_normal, seg1_unit_normal, *current_geom);
                    }
                    else
                    {
                        seg1_unit_normal = offsetter::unit_normal(v1,first_vertex);
                        offset.offset_corner(v1, seg0_unit_normal, seg1_unit_normal, *current_geom);
                        seg0_unit_normal = offsetter::unit_normal(first_vertex, second_vertex);
                        offset.offset_corner(first_vertex, seg1_unit_normal, seg0_unit_normal, *current_geom);
                    }
                }
            }
            else
            {
                // Invalid command
                return false;
            }
        }
        return true;   
    }

    void sort() 
    {
       std::sort(vertices_.begin(), vertices_.end());
    }

    void clear()
    {
        vertices_.clear();
        geometries_.clear();
    }

private:

    std::size_t create_null()
    {
        std::size_t id = geometries_.size();
        geometries_.push_back(clist_null());
        return id;
    }

    void geometry_is_point(clist_vertex_ptr vertex, std::size_t id)
    {
        geometries_[id] = clist_geometry<value_type>(vertex, CLIST_POINT_TYPE);
    }
    
    void geometry_is_point(container_ptr ctr)
    {
        geometry_is_point(ctr->begin(), ctr->get_id());
    }

    void geometry_is_line(clist_vertex_ptr vertex, std::size_t id)
    {
        geometries_[id] = clist_geometry<value_type>(vertex, CLIST_LINE_TYPE);
    }
    
    void geometry_is_line(container_ptr ctr)
    {
        geometry_is_line(ctr->begin(), ctr->get_id());
    }
    
    void geometry_is_ring(clist_vertex_ptr vertex, std::size_t id)
    {
        geometries_[id] = clist_geometry<value_type>(vertex, CLIST_RING_TYPE);
    }
    
    void geometry_is_ring(container_ptr ctr)
    {
        geometry_is_ring(ctr->begin(), ctr->get_id());
    }

    container_ptr new_geometry()
    {
        std::size_t id = create_null();  
        return std::make_shared<container_t>(vertices_, id);
    }
    
    clist_vertex_ptr add_vertex(value_type x, value_type y, std::size_t id, clist_vertex_ptr vertex)
    {
        vertices_.emplace_back(x, y, id);
        clist_vertex_ptr next = &(vertices_.back());
        clist_algo::link_after(vertex, next);
        return next;
    }
    
    clist_vertex_ptr add_before_vertex(value_type x, value_type y, std::size_t id, clist_vertex_ptr vertex)
    {
        vertices_.emplace_back(x, y, id);
        clist_vertex_ptr previous = &(vertices_.back());
        clist_algo::link_before(vertex, previous);
        return previous;
    }
    
    clist_vertex_ptr add_first_vertex(value_type x, value_type y, std::size_t id)
    {
        vertices_.emplace_back(x, y, id);
        clist_vertex_ptr first = &(vertices_.back());
        clist_algo::init_header(first);
        return first;
    }
    
    clist_vertex_ptr remove_last_vertex()
    {
        clist_vertex_ptr last = &(vertices_.back());
        clist_algo::unlink(last);
        vertices_.pop_back();
        return &(vertices_.back());
    }
};

} // end ns geometry

} // end ns mapnik

#endif // MAPNIK_GEOMETRY_CIRCULAR_LIST_HPP
