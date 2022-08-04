#include "casm/configuration/clusterography/impact_neighborhood.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace clust {

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
void add_to_flower_neighborhood(IntegralCluster const &phenomenal,
                                std::set<xtal::UnitCellCoord> &neighborhood,
                                IntegralCluster const &cluster) {
  for (auto const &site : cluster) {
    for (auto const &phenom_site : phenomenal) {
      if (site.sublattice() == phenom_site.sublattice()) {
        xtal::UnitCell trans = phenom_site.unitcell() - site.unitcell();
        for (auto const &tsite : cluster) {
          neighborhood.insert(tsite + trans);
        }
      }
    }
  }
}

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
void add_to_flower_neighborhood(IntegralCluster const &phenomenal,
                                std::set<xtal::UnitCellCoord> &neighborhood,
                                std::set<IntegralCluster> const &orbit) {
  for (auto const &cluster : orbit) {
    add_to_flower_neighborhood(phenomenal, neighborhood, cluster);
  }
}

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
///
/// Notes:
/// - Populate a set of xtal::UnitCellCoord which, if their DoF
///   change, might affect a property associated with the phenomenal
///   cluster.
/// - `orbits` should be the prim periodic orbits associated with
///   non-zero eci.
void add_to_flower_neighborhood(
    IntegralCluster const &phenomenal,
    std::set<xtal::UnitCellCoord> &neighborhood,
    std::vector<std::set<IntegralCluster>> const &orbits) {
  for (auto const &orbit : orbits) {
    add_to_flower_neighborhood(phenomenal, neighborhood, orbit);
  }
}

/// \brief Add cluster sites to neighborhood
void add_to_local_neighborhood(std::set<xtal::UnitCellCoord> &neighborhood,
                               IntegralCluster const &cluster) {
  for (auto const &site : cluster) {
    neighborhood.insert(site);
  }
}

/// \brief Add cluster sites to neighborhood,
///     for each cluster in orbit
void add_to_local_neighborhood(std::set<xtal::UnitCellCoord> &neighborhood,
                               std::set<IntegralCluster> const &orbit) {
  for (auto const &cluster : orbit) {
    add_to_local_neighborhood(neighborhood, cluster);
  }
}

/// \brief Add cluster sites to neighborhood, for each orbit
///
/// Notes:
/// - Populate a set of xtal::UnitCellCoord which, if their DoF
///   change, might affect a property associated with the phenomenal
///   cluster.
/// - `orbits` should be the local-cluster orbits of the
///   phenomenal cluster associated with non-zero eci.
void add_to_local_neighborhood(
    std::set<xtal::UnitCellCoord> &neighborhood,
    std::vector<std::set<IntegralCluster>> const &orbits) {
  for (auto const &orbit : orbits) {
    add_to_local_neighborhood(neighborhood, orbit);
  }
}

}  // namespace clust
}  // namespace CASM
