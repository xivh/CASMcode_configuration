#include "casm/configuration/clusterography/occ_counter.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace clust {

/// \brief Counter over cluster occupations
Counter<std::vector<int>> make_occ_counter(IntegralCluster const &cluster,
                                           xtal::BasicStructure const &prim) {
  std::vector<int> max_occupant_index;
  for (auto const &site : cluster) {
    Index b = site.sublattice();
    max_occupant_index.push_back(prim.basis()[b].occupant_dof().size() - 1);
  }
  return Counter<std::vector<int>>(std::vector<int>(cluster.size(), 0),
                                   max_occupant_index,
                                   std::vector<int>(cluster.size(), 1));
}

}  // namespace clust
}  // namespace CASM
