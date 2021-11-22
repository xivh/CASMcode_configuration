#ifndef CASM_config_SupercellSymInfo
#define CASM_config_SupercellSymInfo

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing application of symmetry in a supercell
struct SupercellSymInfo {
  /// \brief Constructor
  SupercellSymInfo(
      std::shared_ptr<Prim const> const &prim, Superlattice const &superlattice,
      xtal::UnitCellIndexConverter const &unitcell_index_converter,
      xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter);

  /// \brief The subgroup of the prim factor group that leaves
  /// the supercell lattice vectors invariant
  std::shared_ptr<SymGroup const> factor_group;

  /// \brief Describes how sites permute due to translations within the
  /// supercell
  ///
  /// The number of translations is equal the supercell volume (as an integer
  /// multiple of the prim unit cell)
  std::vector<Permutation> translation_permutations;

  /// \brief Describes how sites permute due to supercell factor group
  /// operations.
  ///
  /// There is one element for each element in the supercell factor group.
  std::vector<Permutation> factor_group_permutations;
};

}  // namespace config
}  // namespace CASM

#endif
