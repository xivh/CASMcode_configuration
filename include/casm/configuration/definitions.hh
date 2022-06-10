// The `casm/configuration` module enables apply symmetry to ConfigDoFValues.
//
// Primarily, this purpose of this module is to provide:
//
// - `struct Configuration`: Data structure encapsulating configuration degree
// of freedom (DoF) values and all information necessary to apply symmetry
// operations. This data structure contains:
//   - `ConfigDoFValues dof_values`: The raw DoF values
//   - `std::shared_ptr<Supercell const> supercell`: Specifies all the
//   structural and symmetry information common for all configurations with the
//   same supercell. All members are const.
// - `class SupercellOpIterator`: Class that enables iterating over all symmetry
// operations that are compatible with a particular supercell and applying them
// to a Configuration.
// - Methods for applying symmetry operations:
//   - `ConfigDoFValues &apply(
//          SupercellOpIterator const &op,
//          ConfigDoFValues &dof_values)`;
//   - `ConfigDoFValues copy_apply(
//          SupercellOpIterator const &op,
//          ConfigDoFValues dof_values);`
// - Methods for comparing configurations and finding canonical forms:
//  - `class ConfigCompare`: Provides "less than" comparison of Configuration
//  - `class ConfigIsEquivalent`: Provides "equal to" comparison of
//  Configuration
//  - `is_canonical`, `canonical_form`, `to_canonical`, `from_canonical`:
//     Methods that provide checking and finding Configuration canonical forms
//     (the "greatest" of all symmetrically equivalent configuration).
//
// Allowed dependencies:
// - casm/crystallography
// - casm/clexulator

#ifndef CASM_config_definitions
#define CASM_config_definitions

#include <memory>
#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
}

namespace xtal {
class BasicStructure;
class Lattice;
class SimpleStructure;
class Superlattice;
struct SymOp;
class UnitCell;
class UnitCellCoord;
class UnitCellIndexConverter;
class UnitCellCoordIndexConverter;
}  // namespace xtal

namespace config {
namespace basic_symmetry {
template <typename ElementType>
struct Group;
}

using clexulator::ConfigDoFValues;
using xtal::BasicStructure;
using xtal::Lattice;
using xtal::SimpleStructure;
using xtal::Superlattice;
using xtal::SymOp;
using xtal::UnitCell;
using xtal::UnitCellCoord;

struct Configuration;
struct Prim;
struct PrimSymInfo;
struct Supercell;
struct SupercellSymInfo;
struct UnitCellCoordRep;

typedef long Index;
typedef std::string DoFKey;

/// \brief A basic_symmetry::Group of xtal::SymGroup
typedef basic_symmetry::Group<SymOp> SymGroup;

/// \brief Describes how sublattices are permuted
typedef std::vector<UnitCellCoordRep> BasisPermutationSymGroupRep;

/// \brief One matrix per master group operation
typedef std::vector<Eigen::MatrixXd> GlobalDoFSymGroupRep;

/// \brief One matrix per sublattice
typedef std::vector<Eigen::MatrixXd> LocalDoFSymOpRep;

/// \brief One LocalDoFSymOpRep per master group operation
typedef std::vector<LocalDoFSymOpRep> LocalDoFSymGroupRep;

/// \brief Describes a permutation
///
/// The following convention is used:
///
///        after[i] = before[permutation[i]];
///
typedef std::vector<Index> Permutation;

/// \brief One permutation per sublattice
typedef std::vector<Permutation> OccSymOpRep;

/// \brief One OccSymOpRep per master group operation
typedef std::vector<OccSymOpRep> OccSymGroupRep;

/// \brief Permute container
///
/// Container must support operator[] indexing and copy construction
template <typename Container>
Container &apply(Permutation const &perm, Container &before);

/// \brief Generate permuted copy of a container
///
/// Uses operator[] indexing and copy construction
template <typename Container>
Container copy_apply(Permutation const &perm, Container const &before);

/// \brief Make the inverse permutation
Permutation inverse(Permutation const &perm);

/// \brief Return permutation that is equivalent to applying two permutations
/// sequentially (applied in the order 'first' then 'second')
Permutation combined_permute(Permutation const &first,
                             Permutation const &second);

// --- Inline definitions ---

/// \brief Permute container
///
/// Container must support operator[] indexing and copy construction
template <typename Container>
Container &apply(Permutation const &perm, Container &before) {
  before = copy_apply(perm, before);
  return before;
}

/// \brief Generate permuted copy of a container
///
/// Container must support operator[] indexing and copy construction
template <typename Container>
Container copy_apply(Permutation const &perm, Container const &before) {
  if (before.size() != perm.size()) {
    throw std::runtime_error(
        "Error in copy_apply(Permutation const &, Container const &): "
        "permutation size does not match container size");
  }

  Container after(before);
  for (Index i = 0; i < perm.size(); i++) {
    after[i] = before[perm[i]];
  }
  return after;
}

/// \brief Make the inverse permutation
inline Permutation inverse(Permutation const &perm) {
  Permutation _inverse_perm(perm.size(), 0);
  for (Index i = 0; i < perm.size(); i++) {
    _inverse_perm[perm[i]] = i;
  }
  return _inverse_perm;
}

/// \brief Return permutation that is equivalent to applying two permutations
/// sequentially (applied in the order 'first' then 'second')
inline Permutation combined_permute(Permutation const &first,
                                    Permutation const &second) {
  return copy_apply(second, first);
}

}  // namespace config
}  // namespace CASM

#endif
