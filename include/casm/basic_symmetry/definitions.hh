// The `casm/basic_symmetry` module supports groups, orbits, and algorithms
//
// Primarily, this purpose of this module is to provide:
// - Group: a simple group data structure
//   - Includes multiplication_table, inverse, and group-subgroup relationships
// - [todo] Orbit: an orbit data structure
//   - Includes orbit finding, generating group, equivalence map, sorting
// - [todo] LargeOrbit: a orbit data structure for large orbits
//   - Includes orbit finding, excludes generating group, equivalence map,
//     sorting
//
// Allowed dependencies:
// - CASMcode_global

#ifndef CASM_basic_symmetry_definitions
#define CASM_basic_symmetry_definitions

#include <memory>
#include <vector>

#include "casm/external/Eigen/Dense"

namespace CASM {

namespace basic_symmetry {

typedef long Index;

template <typename ElementType>
struct Group;

/// \brief Specifies the multiplication table for a group
///
/// Defined according to:
///
///     element[k] == element[i] * element[j],
///     k = multiplication_table[i][j]
///
typedef std::vector<std::vector<Index>> MultiplicationTable;

}  // namespace basic_symmetry
}  // namespace CASM

#endif
