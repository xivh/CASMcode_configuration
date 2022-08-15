#ifndef CASM_config_enum_perturbations
#define CASM_config_enum_perturbations

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Make the distinct clusters of sites, taking into account the
///     background configuration symmetry
std::set<std::set<Index>> make_distinct_cluster_sites(
    Configuration const &background,
    std::vector<std::set<std::set<Index>>> const &orbits_as_indices);

/// \brief Make configurations that are distinct occupation perturbations
std::set<Configuration> make_distinct_perturbations(
    Configuration const &background,
    std::set<std::set<Index>> const &distinct_cluster_sites);

/// \brief Make the distinct clusters of sites, taking into account the
///     event group, supercell, and background configuration symmetry
std::set<std::set<Index>> make_distinct_local_cluster_sites(
    Configuration const &background, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group,
    std::vector<std::set<std::set<Index>>> const &local_orbits_as_indices);

/// \brief Make configurations that are distinct local occupation perturbations
std::set<Configuration> make_distinct_local_perturbations(
    Configuration const &background, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group,
    std::set<std::set<Index>> const &distinct_local_cluster_sites);

}  // namespace config
}  // namespace CASM

#endif
