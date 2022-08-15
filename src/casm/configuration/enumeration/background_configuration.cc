#include "casm/configuration/enumeration/background_configuration.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/copy_configuration.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/io/SymInfo_json_io.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"

namespace CASM {
namespace config {

/// \brief Apply occupation to sites
Configuration &apply_occ(Configuration &configuration,
                         std::vector<Index> const &sites,
                         std::vector<int> const &occ) {
  for (int i = 0; i < sites.size(); ++i) {
    configuration.dof_values.occupation[sites[i]] = occ[i];
  }
  return configuration;
}

/// \brief Copy configuration and apply occupation to sites
Configuration copy_apply_occ(Configuration configuration,
                             std::vector<Index> const &sites,
                             std::vector<int> const &occ) {
  apply_occ(configuration, sites, occ);
  return configuration;
}

/// \brief Make the canonical form for a configuration in the context
///    of an occupation event
///
/// \param configuration A configuration
/// \param event_sites Linear sites indices of the cluster of sites that
///     change during the event
/// \param occ_init Initial occupation on sites
/// \param occ_final Final occupation on sites
/// \param begin, end Iterators over the SupercellSymOp consistent with both
///     the supercell of configuration and a local subgroup of the prim factor
///     group (for example a cluster group).
///
/// \return The canonical configuration, allowing either the initial or
///     final occupation on the event sites.
Configuration make_canonical_form(
    Configuration const &configuration, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group) {
  auto begin = std::begin(event_group);
  auto end = std::end(event_group);
  Configuration config_init =
      copy_apply_occ(configuration, event_sites, occ_init);
  Configuration canonical_config_init =
      make_canonical_form(config_init, begin, end);

  Configuration config_final =
      copy_apply_occ(configuration, event_sites, occ_final);
  Configuration canonical_config_final =
      make_canonical_form(config_final, begin, end);

  if (canonical_config_final > canonical_config_init) {
    return canonical_config_final;
  }
  return canonical_config_init;
}

/// \brief Make all configurations, equivalent as infinite crystals under prim
///     factor group operations, that fit into the same supercell, but are
///     inequivalent under the action of a local group.
///
/// \param motif The motif for the background configurations
/// \param supercell The supercell in which to generate distinct background
///     configurations
/// \param event_sites Linear sites indices of the cluster of sites that
///     change during the event
/// \param occ_init Initial occupation on event_sites
/// \param occ_final Final occupation on event_sites
/// \param event_group The SupercellSymOp consistent with both
///     the supercell of configuration and a local subgroup of the prim factor
///     group (for example a cluster group).
///
/// \param The configuration symmetrically equivalent to the background
///     configuration which form symmetrically distinct backgrounds for the
///     event.
std::set<Configuration> make_distinct_background_configurations(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell,
    std::vector<Index> const &event_sites, std::vector<int> const &occ_init,
    std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group) {
  std::vector<Configuration> all =
      make_all_super_configurations(motif, supercell);

  std::set<Configuration> distinct;
  for (Configuration const &configuration : all) {
    distinct.emplace(make_canonical_form(configuration, event_sites, occ_init,
                                         occ_final, event_group));
  }
  return distinct;
}

}  // namespace config
}  // namespace CASM
