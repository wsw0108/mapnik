//
// Created by wsw on 10/23/15.
//

#ifndef MAPTALKS_UTILS_HPP
#define MAPTALKS_UTILS_HPP

// stl
#include <string>
#include <algorithm>

// boost
#include <boost/filesystem.hpp>

// jni
#include "java.hpp"

namespace fs = boost::filesystem;

class maptalks_utils {
public:
  static std::vector<std::string> classpath(fs::path const &dir) {
    std::vector<std::string> files;
    typedef std::vector<fs::path> vec;
    vec paths;
    std::copy(fs::directory_iterator(dir), fs::directory_iterator(),
              std::back_inserter(paths));
    for (vec::const_iterator it(paths.begin()); it != paths.end(); ++it) {
      files.push_back(it->string());
    }
    return files;
  }

  static std::string jstring_to_string(JNIEnv *env, jstring jstr) {
    const char *cstr = env->GetStringUTFChars(jstr, NULL);
    std::string str = cstr;
    env->ReleaseStringUTFChars(jstr, cstr);
    return str;
  }

  static void check_java_exception(JNIEnv *env) {
    if (env->ExceptionCheck()) {
      jthrowable ex = env->ExceptionOccurred();
      env->ExceptionClear();
      jclass class_Throwable = env->FindClass("java/lang/Throwable");
      jmethodID method_Throwable_toString =
          env->GetMethodID(class_Throwable, "toString", "()Ljava/lang/String;");
      jstring jstr =
          (jstring)env->CallObjectMethod(ex, method_Throwable_toString);
      std::string str = maptalks_utils::jstring_to_string(env, jstr);
      throw java_exception(str);
    }
  }
};

#endif // MAPTALKS_UTILS_HPP
