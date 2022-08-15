#ifndef CASM_clust_IntegralClusterOrbitGenerator
#define CASM_clust_IntegralClusterOrbitGenerator

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
struct UnitCellCoordRep;
}

namespace clust {

/// \brief Specifies a particular cluster to use to generate an orbit,
///     and whether sub-cluster orbits should also be generated
///
/// \ingroup Clusterography
struct IntegralClusterOrbitGenerator {
  explicit IntegralClusterOrbitGenerator(IntegralCluster const &_prototype,
                                         bool _include_subclusters = true);

  IntegralCluster prototype;
  bool include_subclusters;
};

inline IntegralClusterOrbitGenerator::IntegralClusterOrbitGenerator(
    IntegralCluster const &_prototype, bool _include_subclusters)
    : prototype(_prototype), include_subclusters(_include_subclusters) {}

}  // namespace clust
}  // namespace CASM

#endif
