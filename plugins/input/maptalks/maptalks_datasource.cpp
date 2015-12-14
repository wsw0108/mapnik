// self
#include "maptalks_datasource.hpp"
#include "java.hpp"
#include "maptalks_featureset.hpp"
#include "maptalks_utils.hpp"

#include <algorithm>
#include <climits>
#include <sstream>

#pragma GCC diagnostic push
#include <mapnik/warning_ignore.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/spirit/include/qi.hpp>
#pragma GCC diagnostic pop
#include <boost/filesystem.hpp>

// mapnik
#include <mapnik/boolean.hpp>
#include <mapnik/unicode.hpp>
#include <mapnik/util/utf_conv_win.hpp>
#include <mapnik/feature.hpp>
#include <mapnik/feature_factory.hpp>
#include <mapnik/feature_kv_iterator.hpp>
#include <mapnik/value_types.hpp>
#include <mapnik/box2d.hpp>
#include <mapnik/debug.hpp>
#include <mapnik/proj_transform.hpp>
#include <mapnik/projection.hpp>
#include <mapnik/util/variant.hpp>
#include <mapnik/util/file_io.hpp>
#include <mapnik/util/geometry_to_ds_type.hpp>
#include <mapnik/make_unique.hpp>
#include <mapnik/geometry_adapters.hpp>
#include <mapnik/json/feature_collection_grammar.hpp>
#include <mapnik/json/extract_bounding_box_grammar_impl.hpp>
#include <mapnik/util/fs.hpp>
#include <mapnik/util/spatial_index.hpp>
#include <mapnik/geom_util.hpp>

namespace fs = boost::filesystem;
using mapnik::datasource;
using mapnik::parameters;

DATASOURCE_PLUGIN(maptalks_datasource)

struct attr_value_converter
    : public boost::static_visitor<mapnik::eAttributeType> {
  mapnik::eAttributeType operator()(mapnik::value_integer /*val*/) const {
    return mapnik::Integer;
  }

  mapnik::eAttributeType operator()(double /*val*/) const {
    return mapnik::Double;
  }

  mapnik::eAttributeType operator()(float /*val*/) const {
    return mapnik::Double;
  }

  mapnik::eAttributeType operator()(bool /*val*/) const {
    return mapnik::Boolean;
  }

  mapnik::eAttributeType operator()(std::string const & /*val*/) const {
    return mapnik::String;
  }

  mapnik::eAttributeType
  operator()(mapnik::value_unicode_string const & /*val*/) const {
    return mapnik::String;
  }

  mapnik::eAttributeType operator()(mapnik::value_null const & /*val*/) const {
    return mapnik::String;
  }
};

maptalks_datasource::maptalks_datasource(parameters const &params)
    : datasource(params), type_(datasource::Vector),
      desc_(maptalks_datasource::name(),
            *params.get<std::string>("encoding", "utf-8")),
      dbname_(*params.get<std::string>("dbname", "")),
      layer_(*params.get<std::string>("layer", "")),
      filter_(*params.get<std::string>("filter", "")),
      page_num_(*params.get<mapnik::value_integer>("page_num", 0)),
      page_size_(*params.get<mapnik::value_integer>("page_size", INT_MAX)),
      extent_(), features_(), tree_(nullptr) {
  MAPNIK_LOG_DEBUG(maptalks_datasource) << "create maptalks datasource";
  this->init(params);
}

void maptalks_datasource::init(mapnik::parameters const &params) {
  if (dbname_.empty()) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: missing <dbname> parameter");
  }

  if (layer_.empty()) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: missing <layer> parameter");
  }

  if (filter_.empty()) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: missing <filter> parameter");
  }

  boost::optional<std::string> home = params.get<std::string>("engine_home");
  if (!home)
    throw mapnik::datasource_exception(
        "maptalks_datasource: missing <engine_home> parameter");

  std::ostringstream option;
  option << "-Dengine.home=" << *home;
  java::option(option.str());

  fs::path lib;
  lib /= *home;
  lib /= "lib";
  std::vector<std::string> files1 = maptalks_utils::classpath(lib);
  for (int i = 0; i < files1.size(); i++) {
    java::classpath(files1[i]);
  }

  fs::path common;
  common /= *home;
  common /= "common";
  std::vector<std::string> files2 = maptalks_utils::classpath(common);
  for (int i = 0; i < files2.size(); i++) {
    java::classpath(files2[i]);
  }

  JavaVM *jvm = java::instance().get_jvm();
  JNIEnv *env = NULL;
  jint rv = JNI_ERR;
  rv = jvm->GetEnv((void **)&env, JNI_VERSION_1_6);
  if (rv == JNI_EDETACHED) {
    JavaVMAttachArgs attachArgs;
    attachArgs.version = JNI_VERSION_1_6;
    attachArgs.name = NULL;
    attachArgs.group = NULL;
    rv = jvm->AttachCurrentThread((void **)&env, &attachArgs);
    MAPNIK_LOG_DEBUG(maptalks_datasource) << "AttachCurrentThread: " << rv;
  }

  if (env == NULL) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: can not get java env");
  }

  jclass class_PathUtils = env->FindClass("org/maptalks/gis/modules/PathUtils");
  jmethodID method_PathUtils_getConfDir = env->GetStaticMethodID(
      class_PathUtils, "getConfDir", "()Ljava/lang/String;");
  jstring path = (jstring)env->CallStaticObjectMethod(
      class_PathUtils, method_PathUtils_getConfDir);
  maptalks_utils::check_java_exception(env);

  jclass class_SettingsStore =
      env->FindClass("org/maptalks/gis/modules/settings/SettingsStore");
  jmethodID constructor_SettingsStore =
      env->GetMethodID(class_SettingsStore, "<init>", "(Ljava/lang/String;)V");
  jobject store =
      env->NewObject(class_SettingsStore, constructor_SettingsStore, path);
  jmethodID method_SettingsStore_readDbSettings =
      env->GetMethodID(class_SettingsStore, "readDbSettings",
                       "()Lcom/alibaba/fastjson/JSONObject;");
  jobject settings =
      env->CallObjectMethod(store, method_SettingsStore_readDbSettings);
  maptalks_utils::check_java_exception(env);

  jclass class_DBDescriptorFactory =
      env->FindClass("org/maptalks/gis/dao/service/DBDescriptorFactory");
  jmethodID method_DBDescriptorFactory_load =
      env->GetStaticMethodID(class_DBDescriptorFactory, "load",
                             "(Lcom/alibaba/fastjson/JSONObject;)V");
  env->CallStaticVoidMethod(class_DBDescriptorFactory,
                            method_DBDescriptorFactory_load, settings);
  maptalks_utils::check_java_exception(env);
  MAPNIK_LOG_DEBUG(maptalks_datasource)
      << "DBDescriptorFactory.load() finished";

  jclass class_MapEnvironment =
      env->FindClass("org/maptalks/gis/dao/service/MapEnvironment");
  jmethodID method_MapEnvironment_getMapService = env->GetStaticMethodID(
      class_MapEnvironment, "getMapService",
      "(Ljava/lang/String;)Lorg/maptalks/gis/dao/service/IMapService;");
  jstring dbname = env->NewStringUTF(dbname_.c_str());
  jobject mapService = env->CallStaticObjectMethod(
      class_MapEnvironment, method_MapEnvironment_getMapService, dbname);
  maptalks_utils::check_java_exception(env);
  if (mapService == NULL) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: can not get map service for db " + dbname_);
  }
  MAPNIK_LOG_DEBUG(maptalks_datasource)
      << "MapEnvironment.getMapService() finished";

  jclass class_MapService = env->GetObjectClass(mapService);
  jmethodID method_MapService_getLayerService = env->GetMethodID(
      class_MapService, "getLayerService",
      "(Ljava/lang/String;)Lorg/maptalks/gis/dao/service/ILayerService;");
  jstring layer = env->NewStringUTF(layer_.c_str());
  jobject layerService = env->CallObjectMethod(
      mapService, method_MapService_getLayerService, layer);
  maptalks_utils::check_java_exception(env);
  if (layerService == NULL) {
    throw mapnik::datasource_exception(
        "maptalks_datasource: can not get layer service for layer " + layer_);
  }
  MAPNIK_LOG_DEBUG(maptalks_datasource)
      << "MapService.getLayerService() finished";

  jclass class_QueryFilter =
      env->FindClass("org/maptalks/gis/dao/service/QueryFilter");
  jmethodID method_QueryFilter_create = env->GetStaticMethodID(
      class_QueryFilter, "create",
      "(Ljava/lang/String;)Lorg/maptalks/gis/dao/service/QueryFilter;");
  jstring filter = env->NewStringUTF(filter_.c_str());
  jobject queryFilter = env->CallStaticObjectMethod(
      class_QueryFilter, method_QueryFilter_create, filter);
  maptalks_utils::check_java_exception(env);
  MAPNIK_LOG_DEBUG(maptalks_datasource) << "QueryFilter.create() finished";

  jclass class_GeoJSON =
      env->FindClass("org/maptalks/gis/core/geojson/GeoJSON");
  jmethodID method_GeoJSON_toString =
      env->GetMethodID(class_GeoJSON, "toString", "()Ljava/lang/String;");
  jclass class_LayerService = env->GetObjectClass(layerService);
  jmethodID method_LayerService_querySpatial =
      env->GetMethodID(class_LayerService, "querySpatial",
                       "(Lorg/maptalks/gis/dao/service/QueryFilter;II)Lorg/"
                       "maptalks/gis/core/geojson/FeatureCollection;");
  jobject collection =
      env->CallObjectMethod(layerService, method_LayerService_querySpatial,
                            queryFilter, page_num_, page_size_);
  maptalks_utils::check_java_exception(env);

  jstring jstr =
      (jstring)env->CallObjectMethod(collection, method_GeoJSON_toString);
  maptalks_utils::check_java_exception(env);
  const char *str = env->GetStringUTFChars(jstr, NULL);

  std::string inline_string(str);
  if (!inline_string.empty()) {
    char const *start = inline_string.c_str();
    char const *end = start + inline_string.size();
    parse_geojson(start, end);
  }

  // TODO: if exception occurred above, this not called
  rv = jvm->DetachCurrentThread();
  MAPNIK_LOG_DEBUG(maptalks_datasource) << "DetachCurrentThread: " << rv;
}

namespace {
using base_iterator_type = char const *;
const mapnik::transcoder geojson_datasource_static_tr("utf8");
const mapnik::json::feature_collection_grammar<base_iterator_type,
                                               mapnik::feature_impl>
    geojson_datasource_static_fc_grammar(geojson_datasource_static_tr);
const mapnik::json::feature_grammar_callback<base_iterator_type,
                                             mapnik::feature_impl>
    geojson_datasource_static_feature_callback_grammar(
        geojson_datasource_static_tr);
const mapnik::json::feature_grammar<base_iterator_type, mapnik::feature_impl>
    geojson_datasource_static_feature_grammar(geojson_datasource_static_tr);
const mapnik::json::extract_bounding_box_grammar<base_iterator_type>
    geojson_datasource_static_bbox_grammar;
}

template <typename Iterator>
void maptalks_datasource::parse_geojson(Iterator start, Iterator end) {
  using boost::spirit::qi::expectation_failure;
  boost::spirit::standard::space_type space;
  mapnik::context_ptr ctx = std::make_shared<mapnik::context_type>();
  std::size_t start_id = 1;

  mapnik::json::default_feature_callback callback(features_);
  Iterator itr = start;

  try {
    bool result = boost::spirit::qi::phrase_parse(
        itr, end, (geojson_datasource_static_fc_grammar)(
                      boost::phoenix::ref(ctx), boost::phoenix::ref(start_id),
                      boost::phoenix::ref(callback)),
        space);
    if (!result || itr != end) {
      itr = start;
      // try parsing as single Feature or single Geometry JSON
      result = boost::spirit::qi::phrase_parse(
          itr, end, (geojson_datasource_static_feature_callback_grammar)(
                        boost::phoenix::ref(ctx), boost::phoenix::ref(start_id),
                        boost::phoenix::ref(callback)),
          space);
      if (!result || itr != end) {
        throw mapnik::datasource_exception("maptalks_datasource: Failed parse "
                                           "GeoJSON file from in-memory "
                                           "string");
      }
    }
  } catch (expectation_failure<char const *> const &ex) {
    itr = start;
    // try parsing as single Feature or single Geometry JSON
    bool result = boost::spirit::qi::phrase_parse(
        itr, end, (geojson_datasource_static_feature_callback_grammar)(
                      boost::phoenix::ref(ctx), boost::phoenix::ref(start_id),
                      boost::phoenix::ref(callback)),
        space);
    if (!result || itr != end) {
      throw mapnik::datasource_exception("maptalks_datasource: Failed parse "
                                         "GeoJSON file from in-memory string");
    }
  }

  using values_container =
      std::vector<std::pair<box_type, std::pair<std::size_t, std::size_t>>>;
  values_container values;
  values.reserve(features_.size());

  std::size_t geometry_index = 0;
  for (mapnik::feature_ptr const &f : features_) {
    mapnik::box2d<double> box = f->envelope();
    if (box.valid()) {
      if (geometry_index == 0) {
        extent_ = box;
        for (auto const &kv : *f) {
          desc_.add_descriptor(mapnik::attribute_descriptor(
              std::get<0>(kv), mapnik::util::apply_visitor(
                                   attr_value_converter(), std::get<1>(kv))));
        }
      } else {
        extent_.expand_to_include(box);
      }
      values.emplace_back(box, std::make_pair(geometry_index, 0));
    }
    ++geometry_index;
  }
  // packing algorithm
  tree_ = std::make_unique<spatial_index_type>(values);
}

maptalks_datasource::~maptalks_datasource() {
  MAPNIK_LOG_DEBUG(maptalks_datasource) << "destroy maptalks datasource";
}

const char *maptalks_datasource::name() { return "maptalks"; }

mapnik::datasource::datasource_t maptalks_datasource::type() const {
  return type_;
}

mapnik::box2d<double> maptalks_datasource::envelope() const { return extent_; }

boost::optional<mapnik::datasource_geometry_t>
maptalks_datasource::get_geometry_type() const {
  boost::optional<mapnik::datasource_geometry_t> result;
  int multi_type = 0;
  unsigned num_features = features_.size();
  for (unsigned i = 0; i < num_features && i < 5; ++i) {
    mapnik::util::to_ds_type(features_[i]->get_geometry());
    if (result) {
      int type = static_cast<int>(*result);
      if (multi_type > 0 && multi_type != type) {
        result.reset(mapnik::datasource_geometry_t::Collection);
        return result;
      }
      multi_type = type;
    }
  }
  return result;
}

mapnik::layer_descriptor maptalks_datasource::get_descriptor() const {
  return desc_;
}

mapnik::featureset_ptr
maptalks_datasource::features(mapnik::query const &q) const {
  // if the query box intersects our world extent then query for features
  mapnik::box2d<double> const &box = q.get_bbox();
  if (extent_.intersects(box)) {
    maptalks_featureset::array_type index_array;
    if (tree_) {
      tree_->query(boost::geometry::index::intersects(box),
                   std::back_inserter(index_array));
      // sort index array to preserve original feature ordering in GeoJSON
      std::sort(index_array.begin(), index_array.end(),
                [](item_type const &item0, item_type const &item1) {
                  return item0.second.first < item1.second.first;
                });
      return std::make_shared<maptalks_featureset>(features_,
                                                   std::move(index_array));
    }
  }
  // otherwise return an empty featureset pointer
  return mapnik::featureset_ptr();
}

mapnik::featureset_ptr
maptalks_datasource::features_at_point(mapnik::coord2d const &pt,
                                       double tol) const {
  mapnik::box2d<double> query_bbox(pt, pt);
  query_bbox.pad(tol);
  mapnik::query q(query_bbox);
  for (auto const &attr_info : desc_.get_descriptors()) {
    q.add_property_name(attr_info.get_name());
  }
  return features(q);
}
