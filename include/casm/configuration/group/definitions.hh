// The `casm/configuration/group` module supports groups, orbits, and algorithms
//
// Primarily, this purpose of this module is to provide:
// - Group: a simple group data structure
//   - Includes multiplication_table, inverse, and group-subgroup relationships
// - make_orbit: function to construct and sort an orbit of objects
// - make_equivalence_map: function to construct an equivalence map
// - make_cyclic_subgroups: function to make cyclic (small) subgroups
// - make_all_subgroups: combines cyclic subgroups to form all subgroups
// - make_invariant_subgroups: make the invariant subgroup of each element of
//   an orbit
//
// Allowed dependencies:
// - CASMcode_global

#ifndef CASM_group_definitions
#define CASM_group_definitions

#include <memory>
#include <vector>

#include "casm/external/Eigen/Dense"

namespace CASM {
namespace group {

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

}  // namespace group
}  // namespace CASM

#endif
