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

/// \brief Make an orbit of periodic clusters
std::set<IntegralCluster> make_prim_periodic_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep);

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_prim_periodic_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep);

// --- Local-cluster orbits ---

/// \brief Copy cluster and apply symmetry operation transformation
IntegralCluster local_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust);

/// \brief Make an orbit of local clusters
std::set<IntegralCluster> make_local_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep);

/// \brief Make groups that leave cluster orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_local_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &phenomenal_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep);

}  // namespace clust
}  // namespace CASM

#endif
