#ifndef CASM_clust_occ_counter
#define CASM_clust_occ_counter

#include <vector>

#include "casm/container/Counter.hh"

namespace CASM {
namespace xtal {
class BasicStructure;
}

namespace clust {
class IntegralCluster;

/// \brief Counter over cluster occupations
Counter<std::vector<int>> make_occ_counter(IntegralCluster const &cluster,
                                           xtal::BasicStructure const &prim);

}  // namespace clust
}  // namespace CASM

#endif
