// The `casm/configuration/irreps` module supports irreducible space
// decomposition
//
// Primarily, this purpose of this module is to provide:
// - IrrepDecomposition: a class to perform irreducible space decompositions
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_sym_info

#ifndef CASM_irreps_definitions
#define CASM_irreps_definitions

#include <memory>
#include <set>
#include <vector>

#include "casm/configuration/sym_info/definitions.hh"
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

}  // namespace irreps
}  // namespace CASM

#endif
