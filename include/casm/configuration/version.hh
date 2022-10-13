#ifndef CASM_config_version
#define CASM_config_version

#include <string>

namespace CASM {
namespace config {

const std::string &version();  // Returns the version defined by the TXT_VERSION
                               // macro at compile time

}
}  // namespace CASM

extern "C" {

/// \brief Return the libcasm_configuration version number
inline const char *casm_config_version() {
  return CASM::config::version().c_str();
}
}

#endif
