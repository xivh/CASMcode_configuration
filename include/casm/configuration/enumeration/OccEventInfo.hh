#ifndef CASM_config_enum_OccEventInfo
#define CASM_config_enum_OccEventInfo

#include "casm/configuration/Prim.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/definitions.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace config {

/// These data structures help construct all the symgroup reps needed
/// for enumeration by providing easier access to the most commonly used
/// data and methods.
/// They are not meant to handle all use cases or be a dependency for
/// other methods in this library.

struct OccEventPrimInfo {
  OccEventPrimInfo(std::shared_ptr<Prim const> const &_prim,
                   occ_events::OccEvent const &_event);

  std::shared_ptr<Prim const> prim;
  occ_events::OccEvent event;

  /// \brief Symgroup rep of prim->sym_info.factor_group
  ///
  /// Use to apply prim factor group elements to `event`
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;

  /// \brief Subgroup of prim factor group that leaves `event` invariant
  std::shared_ptr<SymGroup const> invariant_group;

  /// \brief Symgroup rep of invariant_group.
  ///
  /// Use to transform xtal::UnitCellCoord while leaving the event unchanged,
  /// such as constructing local-cluster orbits.
  std::vector<xtal::UnitCellCoordRep> invariant_group_unitcellcoord_rep;

  /// \brief Make local-cluster orbits from clusters
  std::vector<std::set<clust::IntegralCluster>> make_local_orbits(
      std::set<clust::IntegralCluster> const &local_clusters) const;
};

struct OccEventSupercellInfo {
  OccEventSupercellInfo(
      std::shared_ptr<OccEventPrimInfo const> const &_event_info,
      std::shared_ptr<Supercell const> const &_supercell);

  std::shared_ptr<OccEventPrimInfo const> event_prim_info;
  std::shared_ptr<Supercell const> supercell;

  /// \brief Linear supercell site indices of event sites
  std::vector<Index> sites;

  /// \brief Initial occupation of event sites
  std::vector<int> occ_init;

  /// \brief Final occupation of event sites
  std::vector<int> occ_final;

  /// \brief Symgroup rep of event_info->invariant_group.
  ///
  /// Use to transform a Configuration while leaving the event unchanged,
  /// such as when identifying symmetrically distinct local environments.
  std::vector<SupercellSymOp> supercellsymop_symgroup_rep;

  /// \brief Make canonical local environment configuration
  Configuration make_canonical_form(Configuration const &configuration) const;

  /// \brief Make all configurations, equivalent as infinite crystals under prim
  ///     factor group operations, that fit into the same supercell, but are
  ///     inequivalent under the action of a local group.
  std::set<Configuration> make_distinct_background_configurations(
      Configuration const &motif) const;

  /// \brief Make configurations that are distinct perturbations of local
  /// clusters
  std::set<Configuration> make_distinct_local_perturbations(
      Configuration const &background,
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits) const;

  /// \brief Generate local-cluster orbits and make configurations that are
  ///     distinct perturbations of local clusters
  std::set<Configuration> make_distinct_local_perturbations(
      Configuration const &background,
      std::set<clust::IntegralCluster> const &local_clusters) const;

  /// \brief Make configurations that are distinct perturbations of local
  /// clusters, using all equivalents of motif that fill the supercell
  std::set<Configuration> make_all_distinct_local_perturbations(
      Configuration const &motif,
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits) const;

  /// \brief Generate local-cluster orbits and make configurations that are
  ///     distinct perturbations of local clusters, using all equivalents
  ///     of motif that fill the supercell
  std::set<Configuration> make_all_distinct_local_perturbations(
      Configuration const &motif,
      std::set<clust::IntegralCluster> const &local_clusters) const;

  /// \brief Generate local-cluster orbits and make configurations that are
  ///     distinct perturbations of sites within a cutoff radius of sites
  ///     in the event
  std::set<Configuration> make_all_distinct_local_perturbations(
      Configuration const &motif, double cutoff_radius) const;
};

}  // namespace config
}  // namespace CASM

#endif
