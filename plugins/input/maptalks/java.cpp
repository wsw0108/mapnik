#include "java.hpp"

#include <iostream>
#include <sstream>

std::vector<std::string> java::classpath_;
std::vector<std::string> java::options_;

void java::classpath(const std::string path) {
  typedef std::vector<std::string>::const_iterator iter_type;
  iter_type begin = java::classpath_.begin();
  iter_type end = java::classpath_.end();
  iter_type iter = std::find(begin, end, path);
  if (iter == end) {
    java::classpath_.push_back(path);
  }
}

void java::option(const std::string option) {
  typedef std::vector<std::string>::const_iterator iter_type;
  iter_type begin = java::options_.begin();
  iter_type end = java::options_.end();
  iter_type iter = std::find(begin, end, option);
  if (iter == end) {
    java::options_.push_back(option);
  }
}

JavaVM *java::get_jvm() { return jvm_; }

JNIEnv *java::get_env() { return env_; }

java::java() : jvm_(NULL), env_(NULL) { this->create_jvm(&jvm_, &env_); }

java::~java() { this->destroy_jvm(&jvm_, &env_); }

void java::create_jvm(JavaVM **jvm, JNIEnv **env) {
  std::vector<std::string> paths = java::classpath_;
  std::vector<std::string> options = java::options_;

  int n = 1 + options.size();
  JavaVMOption *opts = new JavaVMOption[n];

  std::ostringstream classpath;
  classpath << "-Djava.class.path=";
  for (int i = 0; i < paths.size(); ++i) {
    if (i != 0) {
#ifdef WIN32
      classpath << ";";
#else
      classpath << ":";
#endif
    }
    classpath << paths[i];
  }
  opts[0].optionString = strdup(classpath.str().c_str());

  for (int i = 0; i < options.size(); ++i) {
    opts[i + 1].optionString = strdup(options[i].c_str());
  }

  JavaVMInitArgs args;
  args.version = JNI_VERSION_1_6;
  args.ignoreUnrecognized = JNI_FALSE;
  args.options = opts;
  args.nOptions = n;
  JavaVM *vm;
  // TODO: use dlopen/LoadLibraryA to avoid link to jvm.so/jvm.dll
  jint z = JNI_CreateJavaVM(&vm, (void **)env, &args);
  *jvm = vm;

  delete[] opts;
}

void java::destroy_jvm(JavaVM **jvm, JNIEnv **env) {
  if (*jvm != NULL) {
    (*jvm)->DestroyJavaVM();
    *jvm = NULL;
    *env = NULL;
  }
}
