#ifndef CASM_config_PrimMagspinInfo
#define CASM_config_PrimMagspinInfo

#include <optional>

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Holds information about continuous magspin DoF or discrete occupants
/// with
///     magspin properties
struct PrimMagspinInfo {
  explicit PrimMagspinInfo(BasicStructure const &structure);

  bool has_continuous_magspin_dof;
  std::optional<std::string> continuous_magspin_key;
  std::optional<std::string> continuous_magspin_flavor;

  bool has_discrete_atomic_magspin_occupants;
  std::optional<std::string> discrete_atomic_magspin_key;
  std::optional<std::string> discrete_atomic_magspin_flavor;
};

}  // namespace config
}  // namespace CASM

#endif
