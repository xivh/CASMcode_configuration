#include "casm/configuration/Configuration.hh"

#include "casm/clexulator/ConfigDoFValuesTools.hh"

namespace CASM {
namespace config {

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell)
    : supercell(_supercell),
      dof_values(clexulator::make_default_config_dof_values(
          _supercell->prim->basicstructure.basis().size(),  // Index N_sublat
          _supercell->superlattice.size(),                  // Index N_volume
          _supercell->prim->global_dof_info,
          _supercell->prim->local_dof_info)) {}

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell,
                             clexulator::ConfigDoFValues const &_dof_values)
    : supercell(_supercell), dof_values(_dof_values) {}

}  // namespace config
}  // namespace CASM
