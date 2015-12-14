#ifndef MAPTALKS_FEATURESET_HPP
#define MAPTALKS_FEATURESET_HPP

// mapnik
#include <mapnik/feature.hpp>

// stl
#include <vector>
#include <deque>

// self
#include "maptalks_datasource.hpp"

class maptalks_featureset : public mapnik::Featureset {
public:
  typedef std::deque<maptalks_datasource::item_type> array_type;

  // this constructor can have any arguments you need
  maptalks_featureset(std::vector<mapnik::feature_ptr> const &features,
                      array_type &&index_array);

  // desctructor
  virtual ~maptalks_featureset();

  // mandatory: you must expose a next() method, called when rendering
  mapnik::feature_ptr next();

private:
  std::vector<mapnik::feature_ptr> const &features_;
  const array_type index_array_;
  array_type::const_iterator index_itr_;
  array_type::const_iterator index_end_;
};

#endif // MAPTALKS_FEATURESET_HPP
