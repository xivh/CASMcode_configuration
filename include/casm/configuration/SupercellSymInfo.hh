#ifndef CASM_config_SupercellSymInfo
#define CASM_config_SupercellSymInfo

#include "casm/configuration/definitions.hh"
#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing application of symmetry in a supercell
struct SupercellSymInfo {
  /// \brief Constructor
  SupercellSymInfo(
      std::shared_ptr<Prim const> const &prim, Superlattice const &superlattice,
      xtal::UnitCellIndexConverter const &unitcell_index_converter,
      xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
      Index max_n_translation_permutations = 100);

  /// \brief The subgroup of the prim factor group that leaves
  /// the supercell lattice vectors invariant
  std::shared_ptr<SymGroup const> factor_group;

  /// \brief Describes how sites permute due to translations within the
  /// supercell
  ///
  /// The number of translations is equal the supercell volume (as an integer
  /// multiple of the prim unit cell). Not populated for large supercells
  /// (n_unitcells > max_n_translation_permutations).
  std::optional<std::vector<sym_info::Permutation>> translation_permutations;

  /// \brief Describes how sites permute due to supercell factor group
  /// operations.
  ///
  /// There is one element for each element in the supercell factor group.
  std::vector<sym_info::Permutation> factor_group_permutations;
};

/// \brief Construct supercell factor group
SymGroup make_factor_group(std::shared_ptr<Prim const> const &prim,
                           Superlattice const &superlattice);

/// \brief Construct a single supercell translation permutation
sym_info::Permutation make_translation_permutation(
    Index translation_index,
    xtal::UnitCellIndexConverter const &ijk_index_converter,
    xtal::UnitCellCoordIndexConverter const &bijk_index_converter);

/// \brief Construct supercell translation permutations
std::vector<sym_info::Permutation> make_translation_permutations(
    xtal::UnitCellIndexConverter const &ijk_index_converter,
    xtal::UnitCellCoordIndexConverter const &bijk_index_converter);

/// \brief Construct supercell factor group permutations
std::vector<sym_info::Permutation> make_factor_group_permutations(
    std::vector<Index> const &head_group_index,
    sym_info::UnitCellCoordSymGroupRep const &unitcellcoord_symgroup_rep,
    xtal::UnitCellCoordIndexConverter const &bijk_index_converter);

}  // namespace config
}  // namespace CASM

#endif
