
#include "casm/configuration/enumeration/OccEventInfo.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/enumeration/ConfigEnumAllOccupations.hh"
#include "casm/configuration/enumeration/background_configuration.hh"
#include "casm/configuration/enumeration/perturbations.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace config {

OccEventPrimInfo::OccEventPrimInfo(std::shared_ptr<Prim const> const &_prim,
                                   occ_events::OccEvent const &_event)
    : prim(_prim),
      event(_event),
      occevent_symgroup_rep(occ_events::make_occevent_symgroup_rep(
          prim->sym_info.factor_group->element, *prim->basicstructure)),
      invariant_group(occ_events::make_occevent_group(
          event, prim->sym_info.factor_group,
          prim->basicstructure->lattice().lat_column_mat(),
          occevent_symgroup_rep)),
      invariant_group_unitcellcoord_rep(
          sym_info::make_unitcellcoord_symgroup_rep(invariant_group->element,
                                                    *prim->basicstructure)) {}

/// \brief Make local-cluster orbits from clusters
std::vector<std::set<clust::IntegralCluster>>
OccEventPrimInfo::make_local_orbits(
    std::set<clust::IntegralCluster> const &local_clusters) const {
  std::set<clust::IntegralCluster> canonical_elements;
  auto const &rep = invariant_group_unitcellcoord_rep;
  for (auto const &cluster : local_clusters) {
    canonical_elements.insert(group::make_canonical_element(
        cluster, rep.begin(), rep.end(), std::less<clust::IntegralCluster>(),
        clust::local_integral_cluster_copy_apply));
  }
  std::vector<std::set<clust::IntegralCluster>> local_orbits;
  for (auto const &cluster : canonical_elements) {
    local_orbits.push_back(make_local_orbit(cluster, rep));
  }
  return local_orbits;
}

OccEventSupercellInfo::OccEventSupercellInfo(
    std::shared_ptr<OccEventPrimInfo const> const &_event_prim_info,
    std::shared_ptr<Supercell const> const &_supercell)
    : event_prim_info(_event_prim_info),
      supercell(_supercell),
      supercellsymop_symgroup_rep(make_local_supercell_symgroup_rep(
          event_prim_info->invariant_group, supercell)) {
  auto cluster_occupation = make_cluster_occupation(event_prim_info->event);
  sites = to_index_vector(cluster_occupation.first,
                          supercell->unitcellcoord_index_converter);
  occ_init = cluster_occupation.second[0];
  occ_final = cluster_occupation.second[1];
}

/// \brief Make canonical local environment configuration
Configuration OccEventSupercellInfo::make_canonical_form(
    Configuration const &configuration) const {
  return CASM::config::make_canonical_form(
      configuration, sites, occ_init, occ_final, supercellsymop_symgroup_rep);
}

/// \brief Make all configurations, equivalent as infinite crystals under prim
///     factor group operations, that fit into the same supercell, but are
///     inequivalent under the action of a local group.
std::set<Configuration>
OccEventSupercellInfo::make_distinct_background_configurations(
    Configuration const &motif) const {
  return CASM::config::make_distinct_background_configurations(
      motif, supercell, sites, occ_init, occ_final,
      supercellsymop_symgroup_rep);
}

/// \brief Make configurations that are distinct perturbations of local clusters
///
/// \param background The background configuration
/// \param local_orbits Local-cluster orbits, generated without consideration
///     of the background configuration. These orbits are broken based on the
///     background configuration symmetry to find all the distinct local
///     environment perturbations.
std::set<Configuration>
OccEventSupercellInfo::make_distinct_local_perturbations(
    Configuration const &background,
    std::vector<std::set<clust::IntegralCluster>> const &local_orbits) const {
  if (background.supercell != supercell) {
    throw std::runtime_error(
        "Error in OccEventSupercellInfo::make_distinct_local_perturbations: "
        "background supercell does not match this supercell");
  }
  auto local_orbits_as_indices = make_orbits_as_indices(
      local_orbits, supercell->unitcellcoord_index_converter);
  auto distinct_local_cluster_sites = make_distinct_local_cluster_sites(
      background, sites, occ_init, occ_final, supercellsymop_symgroup_rep,
      local_orbits_as_indices);
  return CASM::config::make_distinct_local_perturbations(
      background, sites, occ_init, occ_final, supercellsymop_symgroup_rep,
      distinct_local_cluster_sites);
}

/// \brief Make configurations that are distinct perturbations of local clusters
///
/// \param background The background configuration
/// \param local_clusters Local-clusters, used to generate local-orbits
///     without consideration of the background configuration. These orbits
///     are broken based on the background configuration symmetry to find
///     all the distinct local environment perturbations.
///
std::set<Configuration>
OccEventSupercellInfo::make_distinct_local_perturbations(
    Configuration const &background,
    std::set<clust::IntegralCluster> const &local_clusters) const {
  if (background.supercell != supercell) {
    throw std::runtime_error(
        "Error in OccEventSupercellInfo::make_distinct_local_perturbations: "
        "background supercell does not match this supercell");
  }
  return this->make_distinct_local_perturbations(
      background, event_prim_info->make_local_orbits(local_clusters));
}

/// \brief Make configurations that are distinct perturbations of local
/// clusters, using all equivalents of motif that fill the supercell
///
/// \param motif Used to generate distinct background configuration
/// \param local_orbits Local-cluster orbits, generated without consideration
///     of the background configuration. These orbits are broken based on the
///     background configuration symmetry to find all the distinct local
///     environment perturbations.
std::set<Configuration>
OccEventSupercellInfo::make_all_distinct_local_perturbations(
    Configuration const &motif,
    std::vector<std::set<clust::IntegralCluster>> const &local_orbits) const {
  auto distinct_backgrounds =
      this->make_distinct_background_configurations(motif);
  std::set<Configuration> all;
  for (auto const &background : distinct_backgrounds) {
    auto tmp =
        this->make_distinct_local_perturbations(background, local_orbits);
    for (auto const &c : tmp) {
      all.insert(c);
    }
  }
  return all;
}

/// \brief Generate local-cluster orbits and make configurations that are
///     distinct perturbations of local clusters, using all equivalents
///     of motif that fill the supercell
///
/// \param motif Used to generate distinct background configuration
/// \param local_clusters Local-clusters, used to generate local-orbits
///     without consideration of the background configuration. These orbits
///     are broken based on the background configuration symmetry to find
///     all the distinct local environment perturbations.
std::set<Configuration>
OccEventSupercellInfo::make_all_distinct_local_perturbations(
    Configuration const &motif,
    std::set<clust::IntegralCluster> const &local_clusters) const {
  return this->make_all_distinct_local_perturbations(
      motif, event_prim_info->make_local_orbits(local_clusters));
}

/// \brief Generate local-cluster orbits and make configurations that are
///     distinct perturbations of sites within a cutoff radius of sites
///     in the event
std::set<Configuration>
OccEventSupercellInfo::make_all_distinct_local_perturbations(
    Configuration const &motif, double cutoff_radius) const {
  // get sites using cutoff_radius_neighborhood
  clust::CandidateSitesFunction f = clust::cutoff_radius_neighborhood(
      make_cluster(event_prim_info->event), cutoff_radius);
  auto const &converter = supercell->unitcellcoord_index_converter;
  std::vector<xtal::UnitCellCoord> site_coords =
      f(*event_prim_info->prim->basicstructure, clust::dof_sites_filter());
  std::set<Index> site_indices;
  for (auto const &site : site_coords) {
    site_indices.insert(converter(site));
  }

  // get distinct backgrounds
  auto distinct_backgrounds =
      this->make_distinct_background_configurations(motif);

  // for each background, enumerate local occupations
  std::set<Configuration> distinct_local_perturbations;
  for (auto const &background : distinct_backgrounds) {
    ConfigEnumAllOccupations enumerator(background, site_indices);
    while (enumerator.is_valid()) {
      distinct_local_perturbations.emplace(
          config::make_canonical_form(enumerator.value(), this->sites, occ_init,
                                      occ_final, supercellsymop_symgroup_rep));
      enumerator.advance();
    }
  }
  return distinct_local_perturbations;
}

}  // namespace config
}  // namespace CASM
