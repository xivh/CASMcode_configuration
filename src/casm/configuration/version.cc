#include "casm/configuration/version.hh"

using namespace CASM;
using namespace CASM::config;

#ifndef TXT_VERSION
#define TXT_VERSION "unknown"
#endif

const std::string &CASM::config::version() {
  static const std::string &ver = TXT_VERSION;
  return ver;
};
