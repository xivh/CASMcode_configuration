// The `casm/configuration/irreps` module supports irreducible space
// decomposition
//
// Primarily, this purpose of this module is to provide:
// - IrrepDecomposition: a class to perform irreducible space decompositions
//
// Allowed dependencies:
// - CASMcode_global

#ifndef CASM_irreps_definitions
#define CASM_irreps_definitions

#include <iomanip>
#include <memory>
#include <set>
#include <vector>

#include "casm/casm_io/Log.hh"
#include "casm/container/multivector.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace irreps {

typedef long Index;

typedef std::vector<Eigen::MatrixXd> MatrixRep;
typedef std::set<Index> GroupIndices;
typedef std::set<GroupIndices> GroupIndicesOrbit;
typedef std::set<GroupIndicesOrbit> GroupIndicesOrbitSet;

typedef std::set<Index> SubgroupIndices;
typedef std::set<SubgroupIndices> SubgroupOrbit;
typedef std::set<SubgroupOrbit> SubgroupOrbitSet;
typedef std::vector<std::vector<std::vector<Index>>> SubgroupOrbitVec;

inline Eigen::MatrixXd real_I(Index rows, Index cols) {
  return Eigen::MatrixXd::Identity(rows, cols);
}

inline Eigen::MatrixXd real_Zero(Index rows, Index cols) {
  return Eigen::MatrixXd::Zero(rows, cols);
}

inline Eigen::MatrixXcd complex_I(Index rows, Index cols) {
  return Eigen::MatrixXcd::Identity(rows, cols);
}

inline Eigen::MatrixXcd complex_Zero(Index rows, Index cols) {
  return Eigen::MatrixXcd::Zero(rows, cols);
}

inline void append_time(Log &log, int n_newlines) {
  if (log.print()) {
    log.ostream() << " - Time: " << std::setprecision(6) << log.time_s()
                  << " (s)";
    for (int i = 0; i < n_newlines; ++i) {
      log.ostream() << std::endl;
    }
  }
}

}  // namespace irreps
}  // namespace CASM

#endif
