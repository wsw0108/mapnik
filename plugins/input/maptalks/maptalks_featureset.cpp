// mapnik
#include <mapnik/feature.hpp>

// stl
#include <string>
#include <vector>
#include <deque>

// self
#include "maptalks_featureset.hpp"

maptalks_featureset::maptalks_featureset(
    std::vector<mapnik::feature_ptr> const &features, array_type &&index_array)
    : features_(features), index_array_(std::move(index_array)),
      index_itr_(index_array_.begin()), index_end_(index_array_.end()) {}

maptalks_featureset::~maptalks_featureset() {}

mapnik::feature_ptr maptalks_featureset::next() {
  if (index_itr_ != index_end_) {
    maptalks_datasource::item_type const &item = *index_itr_++;
    std::size_t index = item.second.first;
    if (index < features_.size()) {
      return features_.at(index);
    }
  }
  return mapnik::feature_ptr();
}
