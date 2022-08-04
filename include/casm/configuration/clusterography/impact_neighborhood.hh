#ifndef CASM_clust_impact_neighborhood
#define CASM_clust_impact_neighborhood

#include <set>
#include <vector>

namespace CASM {
namespace xtal {
class UnitCellCoord;
}

namespace clust {
class IntegralCluster;

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
void add_to_flower_neighborhood(IntegralCluster const &phenomenal,
                                std::set<xtal::UnitCellCoord> &neighborhood,
                                IntegralCluster const &cluster);

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
void add_to_flower_neighborhood(IntegralCluster const &phenomenal,
                                std::set<xtal::UnitCellCoord> &neighborhood,
                                std::set<IntegralCluster> const &orbit);

/// \brief Using prim periodic translation symmetry, add all sites
///     that share a cluster with phenomenal sites to neighborhood
void add_to_flower_neighborhood(
    IntegralCluster const &phenomenal,
    std::set<xtal::UnitCellCoord> &neighborhood,
    std::vector<std::set<IntegralCluster>> const &orbits);

/// \brief Add cluster sites to neighborhood
void add_to_local_neighborhood(std::set<xtal::UnitCellCoord> &neighborhood,
                               IntegralCluster const &cluster);

/// \brief Add cluster sites to neighborhood,
///     for each cluster in orbit
void add_to_local_neighborhood(std::set<xtal::UnitCellCoord> &neighborhood,
                               std::set<IntegralCluster> const &orbit);

/// \brief Add cluster sites to neighborhood, for each orbit
void add_to_local_neighborhood(
    std::set<xtal::UnitCellCoord> &neighborhood,
    std::vector<std::set<IntegralCluster>> const &orbits);

}  // namespace clust
}  // namespace CASM

#endif
