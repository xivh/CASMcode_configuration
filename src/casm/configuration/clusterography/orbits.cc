#include "casm/configuration/clusterography/orbits.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/group/subgroups.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace clust {

/// \brief Copy cluster and apply symmetry operation transformation
///
/// \param op, Symmetry operation representation to be applied
/// \param clust, Cluster to transform
///
/// \return cluster, sorted and translated to the origin unit cell
///     after applying the symmetry operation transformation
IntegralCluster prim_periodic_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust) {
  if (!clust.size()) {
    return clust;
  }
  clust.sort();
  xtal::UnitCell pos_init = clust[0].unitcell();
  apply(op, clust);
  clust.sort();
  xtal::UnitCell pos_final = clust[0].unitcell();
  xtal::UnitCell trans = pos_init - pos_final;
  for (xtal::UnitCellCoord &element : clust.elements()) {
    element += trans;
  }
  return clust;
}

/// \brief Find translation that leave cluster sites invariant after
///     transformation, up to a permutation
///
/// \param op, Symmetry operation representation to be applied
/// \param clust, Cluster to transform
///
/// \return translation, such that translation * op * clust is a
///     cluster with the same sites as the original clust, up to
///     a permutation
xtal::UnitCell prim_periodic_integral_cluster_frac_translation(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust) {
  if (!clust.size()) {
    return xtal::UnitCell(0, 0, 0);
  }
  clust.sort();
  xtal::UnitCell pos_init = clust[0].unitcell();
  apply(op, clust);
  clust.sort();
  xtal::UnitCell pos_final = clust[0].unitcell();
  return pos_init - pos_final;
}

/// \brief Make an orbit of periodic clusters
///
/// \param orbit_element One cluster in the orbit
/// \param unitcellcoordrep Symmetry group representation (as
/// xtal::UnitCellCoordRep)
///
std::set<IntegralCluster> make_prim_periodic_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep) {
  return group::make_orbit(orbit_element, unitcellcoordrep.begin(),
                           unitcellcoordrep.end(), std::less<IntegralCluster>(),
                           prim_periodic_integral_cluster_copy_apply);
}

/// \brief Make groups that leave cluster orbit elements invariant
///
/// \param orbit A cluster orbit
/// \param factor_group The factor group used to generate the orbit.
/// \param lat_column_mat The 3x3 matrix whose columns are the lattice vectors
/// \param unitcellcoordrep Symmetry group representation (as
/// xtal::UnitCellCoordRep)
///
/// \returns Cluster invariant groups, where prim_periodic_cluster_groups[i], is
///     the SymGroup whose operations leave the sites of the i-th cluster in the
///     orbit invariant (up to a permutation).
std::vector<std::shared_ptr<SymGroup const>> make_prim_periodic_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep) {
  // The indices eq_map[i] are the indices of the group
  // elements transform the first element in the orbit into the
  // i-th element in the orbit.
  std::vector<std::vector<Index>> eq_map = group::make_equivalence_map(
      orbit, unitcellcoordrep.begin(), unitcellcoordrep.end(),
      prim_periodic_integral_cluster_copy_apply);

  // The indices subgroup_indices[i] are the indices of the group
  // elements which leave orbit element i invariant (up to a translation).
  std::vector<group::SubgroupIndices> subgroup_indices =
      group::make_invariant_subgroups(eq_map, *factor_group);

  // The group cluster_groups[i] contains the SymOp corresponding to
  // subgroup_indices[i] and including the translation which keeps
  // the i-th cluster invariant
  std::vector<std::shared_ptr<SymGroup const>> cluster_groups;
  auto orbit_it = orbit.begin();
  auto subgroup_indices_it = subgroup_indices.begin();
  auto subgroup_indices_end = subgroup_indices.end();

  // return xtal::SymOp, translation * factor_group->element[j], which leaves
  // *orbit_it invariant
  auto make_cluster_group_element = [&](Index j) {
    return xtal::SymOp(
               Eigen::Matrix3d::Identity(),
               lat_column_mat * prim_periodic_integral_cluster_frac_translation(
                                    unitcellcoordrep[j], *orbit_it)
                                    .cast<double>(),
               false) *
           factor_group->element[j];
  };

  while (subgroup_indices_it != subgroup_indices_end) {
    std::vector<xtal::SymOp> cluster_group_elements;
    for (Index j : *subgroup_indices_it) {
      cluster_group_elements.push_back(make_cluster_group_element(j));
    }
    cluster_groups.emplace_back(std::make_shared<SymGroup>(
        factor_group, cluster_group_elements, *subgroup_indices_it));
    ++subgroup_indices_it;
    ++orbit_it;
  }
  return cluster_groups;
}

// --- Local-cluster orbits ---

/// \brief Copy cluster and apply symmetry operation transformation
///
/// \param op, Symmetry operation representation to be applied
/// \param clust, Cluster to transform
///
/// \return cluster, sorted after applying the symmetry operation
///     transformation
IntegralCluster local_integral_cluster_copy_apply(
    xtal::UnitCellCoordRep const &op, IntegralCluster clust) {
  if (!clust.size()) {
    return clust;
  }
  apply(op, clust);
  clust.sort();
  return clust;
}

/// \brief Make an orbit of local clusters
///
/// \param orbit_element One cluster in the orbit
/// \param unitcellcoordrep Symmetry group representation (as
/// xtal::UnitCellCoordRep)
///
std::set<IntegralCluster> make_local_orbit(
    IntegralCluster const &orbit_element,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep) {
  return group::make_orbit(orbit_element, unitcellcoordrep.begin(),
                           unitcellcoordrep.end(), std::less<IntegralCluster>(),
                           local_integral_cluster_copy_apply);
}

/// \brief Make groups that leave cluster orbit elements invariant
///
/// \param orbit A cluster orbit
/// \param head_group The phenomenal cluster group used to generate the orbit.
/// \param unitcellcoordrep Symmetry group representation (as
/// xtal::UnitCellCoordRep)
///
/// \returns Cluster invariant groups, where prim_periodic_cluster_groups[i], is
///     the SymGroup whose operations leave the sites of the i-th cluster in the
///     orbit invariant (up to a permutation).
std::vector<std::shared_ptr<SymGroup const>> make_local_cluster_groups(
    std::set<IntegralCluster> const &orbit,
    std::shared_ptr<SymGroup const> const &phenomenal_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoordrep) {
  // The indices eq_map[i] are the indices of the group
  // elements transform the first element in the orbit into the
  // i-th element in the orbit.
  std::vector<std::vector<Index>> eq_map = group::make_equivalence_map(
      orbit, unitcellcoordrep.begin(), unitcellcoordrep.end(),
      local_integral_cluster_copy_apply);

  // The indices subgroup_indices[i] are the indices of the group
  // elements which leave orbit element i invariant.
  std::vector<group::SubgroupIndices> subgroup_indices =
      make_invariant_subgroups(eq_map, *phenomenal_group);

  // The group cluster_groups[i] contains the SymOp corresponding to
  // subgroup_indices[i].
  std::vector<std::shared_ptr<SymGroup const>> cluster_groups;
  auto orbit_it = orbit.begin();
  auto subgroup_indices_it = subgroup_indices.begin();
  auto subgroup_indices_end = subgroup_indices.end();
  while (subgroup_indices_it != subgroup_indices_end) {
    std::vector<xtal::SymOp> cluster_group_elements;
    for (Index j : *subgroup_indices_it) {
      cluster_group_elements.push_back(phenomenal_group->element[j]);
    }
    cluster_groups.emplace_back(std::make_shared<SymGroup>(
        phenomenal_group, cluster_group_elements, *subgroup_indices_it));
    ++subgroup_indices_it;
    ++orbit_it;
  }
  return cluster_groups;
}

}  // namespace clust
}  // namespace CASM
