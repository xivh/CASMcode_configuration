#ifndef CASM_occ_events_LexicographicalCompare
#define CASM_occ_events_LexicographicalCompare

#include <algorithm>

#include "casm/global/eigen.hh"

namespace CASM {
namespace occ_events {

/// \brief Lexicographically compare Eigen::VectorXi
struct LexicographicalCompare {
  bool operator()(Eigen::VectorXi const &A, Eigen::VectorXi const &B) const {
    return std::lexicographical_compare(A.data(), A.data() + A.size(), B.data(),
                                        B.data() + B.size());
  }
};

}  // namespace occ_events
}  // namespace CASM

#endif
