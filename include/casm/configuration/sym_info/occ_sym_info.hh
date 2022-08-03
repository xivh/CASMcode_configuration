#ifndef CASM_sym_info_occ
#define CASM_sym_info_occ

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace sym_info {

struct OccSymInfo {
  OccSymInfo(std::vector<xtal::SymOp> const &group_elements,
             xtal::BasicStructure const &basicstructure);

  /// \brief True if any occupation DoF
  bool has_occupation_dofs;

  /// \brief True if any permutation in occ_symgroup_rep is non-trivial
  bool has_aniso_occs;

  /// \brief Permutations describe occupant index transformation under symmetry
  ///
  /// Usage:
  /// \code
  /// Index occupant_index_after = atom_position_symgroup_rep
  ///                                  .at(group_element_index)
  ///                                  .at(sublattice_index_before)
  ///                                  .at(occupant_index_before);
  /// \endcode
  ///
  /// Note:
  /// - This describes cases such as discrete molecular orientations or
  /// occupants with spin where a symmetry operation may transform one discrete
  /// occupant into another *before* permutating among sites.
  OccSymGroupRep occ_symgroup_rep;

  /// \brief Permutations describe atom position index transformation under
  /// symmetry
  ///
  /// Usage:
  /// \code
  /// Index atom_position_index_after = atom_position_symgroup_rep
  ///                                       .at(group_element_index)
  ///                                       .at(sublattice_index_before)
  ///                                       .at(occupant_index_before);
  /// \endcode
  ///
  /// Note:
  /// - This describes symmetry operations transforming molecules, resulting in
  /// permutation of atoms among the symmetrically equivalent atom positions in
  /// the molecule.
  AtomPositionSymGroupRep atom_position_symgroup_rep;
};

}  // namespace sym_info
}  // namespace CASM

#endif
