#ifndef CASM_clust_orbits
#define CASM_clust_orbits

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/group/orbits.hh"

namespace CASM {
namespace clust {

/// \brief Make an orbit of periodic clusters
///
/// \param orbit_element One local cluster in the orbit
/// \param group_begin,group_end Iterators over xtal::UnitCellCoordRep
///
template <typename OrbitElementType, typename GroupElementIt,
          typename CopyApplyType>
std::set<IntegralCluster> make_prim_periodic_orbit(
    IntegralCluster const &orbit_element, GroupElementIt group_begin,
    GroupElementIt group_end) {
  auto copy_apply_f = [](xtal::UnitCellCoordRep const &op,
                         IntegralCluster clust) {
    if (!clust.size()) {
      return clust;
    }
    apply(op, clust);
    clust.sort();
    xtal::UnitCell pos = clust[0].unitcell();
    for (UnitCellCoord &element : clust.elements()) {
      element -= pos;
    }
    return clust;
  };
  return make_orbit(orbit_element, group_begin, group_end,
                    std::less<IntegralCluster>, copy_apply_f);
  return orbit;
};

/// \brief Make an orbit of local clusters
///
/// \param orbit_element One local cluster in the orbit
/// \param group_begin,group_end Iterators over a group of
/// xtal::UnitCellCoordRep
///
template <typename OrbitElementType, typename GroupElementIt,
          typename CopyApplyType>
std::set<IntegralCluster> make_local_orbit(IntegralCluster const &orbit_element,
                                           GroupElementIt group_begin,
                                           GroupElementIt group_end) {
  auto copy_apply_f = [](xtal::UnitCellCoordRep const &op,
                         IntegralCluster clust) {
    if (!clust.size()) {
      return clust;
    }
    apply(op, clust);
    clust.sort();
    return clust;
  };
  return make_orbit(orbit_element, group_begin, group_end,
                    std::less<IntegralCluster>, copy_apply_f);
  return orbit;
};

}  // namespace clust
}  // namespace CASM

#endif
