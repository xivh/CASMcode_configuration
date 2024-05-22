#ifndef CASM_clust_orbits
#define CASM_clust_orbits

#include <memory>
#include <set>
#include <vector>

#include "casm/configuration/clusterography/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace clust {

class IntegralCluster;

/// \brief Copy cluster and apply symmetry operation transformation
IntegralCluster prim_periodic_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust);

/// \brief Find translation that leave cluster sites invariant after
///     transformation, up to a permutation
xtal::UnitCell prim_periodic_integral_cluster_frac_translation(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust);

/// \brief Make an orbit of clusters, with periodic symmetry of a prim
std::set<IntegralCluster> make_prim_periodic_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make equivalence map of factor group indices for an orbit of
/// clusters,
///     with periodic symmetry of a prim
std::vector<std::vector<Index>> make_prim_periodic_equivalence_map_indices(
    std::set<IntegralCluster> const &orbit,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make equivalence map for an orbit of clusters, with periodic
///     symmetry of a prim
std::vector<std::vector<xtal::SymOp>> make_prim_periodic_equivalence_map(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &symgroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Return xtal::SymOp that leaves phenomenal invariant, and is a
///     combination of a factor group operation and a lattice translation
xtal::SymOp make_cluster_group_element(
    IntegralCluster const &phenomenal, Eigen::Matrix3d const &lat_column_mat,
    xtal::SymOp const &factor_group_op,
    xtal::UnitCellCoordRep const &unitcellcoord_rep);

/// \brief Find translation necessary to construct an equivalence map operation
xtal::UnitCell equivalence_map_translation(xtal::UnitCellCoordRep const &op,
                                           IntegralCluster prototype,
                                           IntegralCluster equivalent);

/// \brief Return xtal::SymOp that maps prototype to equivalent, and is a
///     combination of a factor group operation and a lattice translation
xtal::SymOp make_equivalence_map_op(
    IntegralCluster const &prototype, IntegralCluster const &equivalent,
    Eigen::Matrix3d const &lat_column_mat, xtal::SymOp const &factor_group_op,
    xtal::UnitCellCoordRep const &unitcellcoord_rep);

/// \brief Determine the translations used to generate known
///     phenomenal clusters from a known prototype and
///     generating ops
std::vector<xtal::UnitCell> make_phenomenal_generating_translations(
    clust::IntegralCluster const &prototype,
    std::vector<clust::IntegralCluster> const &phenomenal_clusters,
    std::vector<Index> const &equivalent_generating_op_indices,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &symggroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make the group which leaves a cluster invariant
std::shared_ptr<SymGroup const> make_cluster_group(
    IntegralCluster cluster, std::shared_ptr<SymGroup const> const &symggroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make orbits of clusters, with periodic symmetry of a prim
std::vector<std::set<IntegralCluster>> make_prim_periodic_orbits(
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
    SiteFilterFunction site_filter, std::vector<double> const &max_length,
    std::vector<IntegralClusterOrbitGenerator> const &custom_generators);

/// \brief Convert orbits of IntegralCluster to orbits of linear site
///     indices in a supercell
std::vector<std::set<std::set<Index>>> make_orbits_as_indices(
    std::vector<std::set<IntegralCluster>> const &orbits,
    xtal::UnitCellCoordIndexConverter const &converter);

// --- Local-cluster orbits ---

/// \brief Minimal "equivalents info" specify the phenomenal
///     clusters of equivalent local basis sets
struct EquivalentsInfo {
  /// \brief The phenomenal clusters
  std::vector<clust::IntegralCluster> phenomenal_clusters;

  /// \brief Indices of the factor group operations that
  ///     generate the equivalent phenomenal clusters
  std::vector<Index> equivalent_generating_op_indices;
};

/// \brief Copy cluster and apply symmetry operation transformation
IntegralCluster local_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust);

/// \brief Make an orbit of local clusters
std::set<IntegralCluster> make_local_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make equivalence map of phenomenal group indices for an orbit of
/// local
///     clusters
std::vector<std::vector<Index>> make_local_equivalence_map_indices(
    std::set<IntegralCluster> const &orbit,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make equivalence map for an orbit of local clusters
std::vector<std::vector<xtal::SymOp>> make_local_equivalence_map(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &phenomenal_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_local_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &phenomenal_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make the group that leaves a local cluster invariant
std::shared_ptr<SymGroup const> make_local_cluster_group(
    IntegralCluster cluster,
    std::shared_ptr<SymGroup const> const &phenomenal_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make local-cluster orbits
std::vector<std::set<IntegralCluster>> make_local_orbits(
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
    SiteFilterFunction site_filter, std::vector<double> const &max_length,
    std::vector<IntegralClusterOrbitGenerator> const &custom_generators,
    IntegralCluster const &phenomenal, std::vector<double> const &cutoff_radius,
    bool include_phenomenal_sites = false);

}  // namespace clust
}  // namespace CASM

#endif
