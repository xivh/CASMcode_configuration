#ifndef CASM_clust_orbits
#define CASM_clust_orbits

#include <set>
#include <vector>

#include "casm/configuration/clusterography/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
struct UnitCellCoordRep;
}

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

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make the group which leaves a cluster invariant
std::shared_ptr<SymGroup const> make_cluster_group(
    IntegralCluster cluster,
    std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make orbits of clusters, with periodic symmetry of a prim
std::vector<std::set<IntegralCluster>> make_prim_periodic_orbits(
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
    SiteFilterFunction site_filter, std::vector<double> const &max_length,
    std::vector<IntegralClusterOrbitGenerator> const &custom_generators);

// --- Local-cluster orbits ---

/// \brief Copy cluster and apply symmetry operation transformation
IntegralCluster local_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust);

/// \brief Make an orbit of local clusters
std::set<IntegralCluster> make_local_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_local_cluster_groups(
    std::set<IntegralCluster> const &orbit,
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
