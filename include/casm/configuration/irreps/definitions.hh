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
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"

namespace CASM {
namespace irreps {

typedef long Index;
typedef std::vector<Eigen::MatrixXd> MatrixRep;
typedef std::set<Index> GroupIndices;
typedef std::set<GroupIndices> GroupIndicesOrbit;
typedef std::set<GroupIndicesOrbit> GroupIndicesOrbitSet;

}  // namespace irreps
}  // namespace CASM

#endif
