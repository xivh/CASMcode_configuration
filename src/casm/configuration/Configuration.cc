#include "casm/configuration/Configuration.hh"

#include "casm/clexulator/ConfigDoFValuesTools.hh"
#include "casm/configuration/ConfigCompare.hh"
#include "casm/configuration/ConfigIsEquivalent.hh"
#include "casm/configuration/SupercellSymOp.hh"

namespace CASM {
namespace config {

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell)
    : supercell(_supercell),
      dof_values(clexulator::make_default_config_dof_values(
          _supercell->prim->basicstructure.basis().size(),
          _supercell->superlattice.size(), _supercell->prim->global_dof_info,
          _supercell->prim->local_dof_info)) {}

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell,
                             clexulator::ConfigDoFValues const &_dof_values)
    : supercell(_supercell), dof_values(_dof_values) {}

/// \brief Less than comparison of Configuration
///
/// - Must have the same Supercell
bool Configuration::operator<(Configuration const &rhs) const {
  ConfigCompare compare(*this);
  return compare(rhs);
}

/// \brief Equality comparison of Configuration
///
/// - Must have the same Supercell
/// - Checks that all DoF are the same, within tolerance
bool Configuration::eq_impl(Configuration const &rhs) const {
  if (supercell != rhs.supercell) {
    return false;
  }
  double xtal_tol = supercell->prim->basicstructure.lattice().tol();
  ConfigIsEquivalent equal_to(*this, xtal_tol);
  return equal_to(rhs);
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration &apply(SupercellSymOp const &op, Configuration &configuration) {
  apply(op, configuration.dof_values);
  return configuration;
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration copy_apply(SupercellSymOp const &op,
                         Configuration configuration) {
  return apply(op, configuration);
}

}  // namespace config
}  // namespace CASM
