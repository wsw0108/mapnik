#ifndef MAPTALKS_JAVA_HPP
#define MAPTALKS_JAVA_HPP

// mapnik
#include <mapnik/util/singleton.hpp>
#include <mapnik/util/noncopyable.hpp>

// stl
#include <string>
#include <vector>

// jni
#include <jni.h>

class java_exception : public std::exception {
public:
  java_exception(std::string const &message) : message_(message) {}

  ~java_exception() throw() {}

  virtual const char *what() const throw() { return message_.c_str(); }

private:
  std::string message_;
};

class java : public mapnik::singleton<java, mapnik::CreateStatic>,
             private mapnik::util::noncopyable {
  friend class mapnik::CreateStatic<java>;

public:
  static void classpath(const std::string path);
  static void option(const std::string option);
  JavaVM *get_jvm();
  JNIEnv *get_env();

private:
  static std::vector<std::string> classpath_;
  static std::vector<std::string> options_;

private:
  java();
  ~java();
  void create_jvm(JavaVM **jvm, JNIEnv **env);
  void destroy_jvm(JavaVM **jvm, JNIEnv **env);
  JavaVM *jvm_;
  JNIEnv *env_;
};

#endif // MAPTALKS_JAVA_HPP
