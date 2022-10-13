// The `casm/configuration/sym_info` module enables generating
// symmetry representations.
//
// Primarily, this purpose of this module is to provide:
//
// - Symmetry representations
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_crystallography
// - CASMcode_configuration/group

#ifndef CASM_sym_info_definitions
#define CASM_sym_info_definitions

#include <map>
#include <memory>
#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace group {
template <typename ElementType>
struct Group;
}

namespace xtal {
class BasicStructure;
class Lattice;
struct SymOp;
class UnitCell;
class UnitCellCoord;
struct UnitCellCoordRep;
}  // namespace xtal

namespace sym_info {

using xtal::BasicStructure;
using xtal::Lattice;
using xtal::SymOp;
using xtal::UnitCell;
using xtal::UnitCellCoord;
using xtal::UnitCellCoordRep;

typedef long Index;
typedef std::string DoFKey;

/// \brief A group::Group of xtal::SymOp
typedef group::Group<SymOp> SymGroup;

/// \brief Describes how sites are permuted
typedef std::vector<UnitCellCoordRep> UnitCellCoordSymGroupRep;

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

/// \brief One permutation per sublattice, per occupant
typedef std::vector<std::vector<Permutation>> AtomPositionSymOpRep;

/// \brief One AtomPositionSymOpRep per master group operation
typedef std::vector<AtomPositionSymOpRep> AtomPositionSymGroupRep;

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

/// \brief Return a permutation matrix, P, such that `P * before`
///      is equivalent to `apply(perm, before)`
Eigen::MatrixXd as_matrix(Permutation const &perm);

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

/// \brief Return a permutation matrix, P, such that `P * before`
///      is equivalent to `apply(perm, before)`
inline Eigen::MatrixXd as_matrix(Permutation const &perm) {
  Eigen::MatrixXd P(perm.size(), perm.size());
  P.setZero();
  for (Index i = 0; i < perm.size(); ++i) {
    P(i, perm[i]) = 1.0;
  }
  return P;
}

}  // namespace sym_info
}  // namespace CASM

#endif
