#ifndef CASM_config_dof_space_functions
#define CASM_config_dof_space_functions

#include <map>

#include "casm/clexulator/DoFSpace.hh"

namespace CASM {
namespace clexulator {
struct DoFSpace;
}
namespace config {

/// Removes specified occupation modes from the DoFSpace basis, by sublattice
clexulator::DoFSpace exclude_default_occ_modes_by_sublattice(
    clexulator::DoFSpace const &dof_space,
    std::map<int, int> sublattice_index_to_default_occ);

/// Removes specified occupation modes from the DoFSpace basis, by supercell
/// site index
clexulator::DoFSpace exclude_default_occ_modes_by_site(
    clexulator::DoFSpace const &dof_space,
    std::map<Index, int> site_index_to_default_occ);

/// Removes specified occupation modes from the DoFSpace basis
clexulator::DoFSpace exclude_default_occ_modes(
    clexulator::DoFSpace const &dof_space_in,
    bool include_default_occ_modes = false,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ =
        std::nullopt,
    std::optional<std::map<Index, int>> site_index_to_default_occ =
        std::nullopt);

/// Removes homogeneous modes from the DoFSpace basis
clexulator::DoFSpace exclude_homogeneous_mode_space(
    clexulator::DoFSpace const &dof_space_in,
    std::optional<bool> exclude_homogeneous_modes = std::nullopt);

}  // namespace config
}  // namespace CASM

#endif
