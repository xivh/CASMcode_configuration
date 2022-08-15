#ifndef CASM_config_enum_background_configuration
#define CASM_config_enum_background_configuration

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// Notes:
///
/// Given an event, background configuration motif, and supercell,
/// want to generate unique perturbations of the local environment
/// of the event.
///
/// Steps:
/// 1) Generate distinct background configurations. These are the
///    configurations that are symmetrically equal to the background
///    configuration w.r.t. the prim factor group in the infinite crystal,
///    but distinct given the position of the event in a particular
///    supercell.
/// 2) Generate distinct local-cluster sites, as indices in the supercell,
///    given the event, background configuration, and local-cluster orbits
///    generated for the event w.r.t. the prim factor group.
/// 3) Generate configurations which represent perturbations of the local
///    environment, and make canonical w.r.t. to the event in order to
///    keep only the distinct perturbations.
///

/// \brief Apply occupation to sites
Configuration &apply_occ(Configuration &configuration,
                         std::vector<Index> const &sites,
                         std::vector<int> const &occ);

/// \brief Copy configuration and apply occupation to sites
Configuration copy_apply_occ(Configuration configuration,
                             std::vector<Index> const &sites,
                             std::vector<int> const &occ);

/// \brief Make the canonical form for a configuration in the context
///    of an occupation event
Configuration make_canonical_form(
    Configuration const &configuration, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group);

/// \brief Make all configurations, equivalent as infinite crystals under prim
///     factor group operations, that fit into the same supercell, but are
///     inequivalent under the action of a local group.
std::set<Configuration> make_distinct_background_configurations(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell,
    std::vector<Index> const &event_sites, std::vector<int> const &occ_init,
    std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group);

}  // namespace config
}  // namespace CASM

#endif
