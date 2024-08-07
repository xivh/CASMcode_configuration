#ifndef CASM_config_ConfigCompare
#define CASM_config_ConfigCompare

#include <utility>

#include "casm/configuration/ConfigIsEquivalent.hh"

namespace CASM {
namespace config {

/// \brief Class for less than comparison of Configurations
class ConfigCompare {
 public:
  explicit ConfigCompare(ConfigIsEquivalent const &_eq) : m_eq(_eq) {}
  explicit ConfigCompare(Configuration const &_config,
                         std::set<std::string> const &_which_dofs)
      : m_eq(_config, _which_dofs) {}

  template <typename... Args>
  bool operator()(Args &&...args) const {
    if (m_eq(std::forward<Args>(args)...)) {
      return false;
    }
    return m_eq.is_less();
  }

 private:
  ConfigIsEquivalent m_eq;
};

}  // namespace config
}  // namespace CASM

#endif
