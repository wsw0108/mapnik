#ifndef MAPTALKS_DATASOURCE_HPP
#define MAPTALKS_DATASOURCE_HPP

// mapnik
#include <mapnik/datasource.hpp>
#include <mapnik/params.hpp>
#include <mapnik/query.hpp>
#include <mapnik/feature.hpp>
#include <mapnik/box2d.hpp>
#include <mapnik/coord.hpp>
#include <mapnik/feature_layer_desc.hpp>
#include <mapnik/unicode.hpp>

#pragma GCC diagnostic push
#include <mapnik/warning_ignore.hpp>
#include <boost/optional.hpp>
#include <boost/version.hpp>
#include <boost/geometry/index/rtree.hpp>
#pragma GCC diagnostic pop

// stl
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <deque>

template <std::size_t Max, std::size_t Min>
struct geojson_linear : boost::geometry::index::linear<Max, Min> {};

namespace boost {
namespace geometry {
namespace index {
namespace detail {
namespace rtree {

template <std::size_t Max, std::size_t Min>
struct options_type<geojson_linear<Max, Min>> {
  using type =
      options<geojson_linear<Max, Min>, insert_default_tag,
              choose_by_content_diff_tag, split_default_tag, linear_tag,
#if BOOST_VERSION >= 105700
              node_variant_static_tag>;
#else
              node_s_mem_static_tag>;

#endif
};
}
}
}
}
}

class maptalks_datasource : public mapnik::datasource {
public:
  using box_type = mapnik::box2d<double>;
  using item_type = std::pair<box_type, std::pair<std::size_t, std::size_t>>;
  using spatial_index_type =
      boost::geometry::index::rtree<item_type, geojson_linear<16, 4>>;

  // constructor
  // arguments must not change
  maptalks_datasource(mapnik::parameters const &params);

  // destructor
  virtual ~maptalks_datasource();

  // mandatory: type of the plugin, used to match at runtime
  mapnik::datasource::datasource_t type() const;

  // mandatory: name of the plugin
  static const char *name();

  // mandatory: function to query features by box2d
  // this is called when rendering, specifically in feature_style_processor.hpp
  mapnik::featureset_ptr features(mapnik::query const &q) const;

  // mandatory: function to query features by point (coord2d)
  // not used by rendering, but available to calling applications
  mapnik::featureset_ptr features_at_point(mapnik::coord2d const &pt,
                                           double tol = 0) const;

  // mandatory: return the box2d of the datasource
  // called during rendering to determine if the layer should be processed
  mapnik::box2d<double> envelope() const;

  // mandatory: optionally return the overal geometry type of the datasource
  boost::optional<mapnik::datasource_geometry_t> get_geometry_type() const;

  // mandatory: return the layer descriptor
  mapnik::layer_descriptor get_descriptor() const;

private:
  void init(mapnik::parameters const &params);

  template <typename Iterator> void parse_geojson(Iterator start, Iterator end);

  mapnik::datasource::datasource_t type_;
  mapnik::layer_descriptor desc_;
  std::string dbname_;
  std::string layer_;
  std::string filter_;
  mapnik::value_integer page_num_;
  mapnik::value_integer page_size_;
  mapnik::box2d<double> extent_;
  std::vector<mapnik::feature_ptr> features_;
  std::unique_ptr<spatial_index_type> tree_;
};

#endif // MAPTALKS_DATASOURCE_HPP
