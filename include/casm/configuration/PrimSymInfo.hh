#ifndef CASM_config_PrimSymInfo
#define CASM_config_PrimSymInfo

#include "casm/configuration/definitions.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/sym_info/definitions.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing how prim DoF values transform under
/// application of symmetry
struct PrimSymInfo {
  /// \brief Constructor
  PrimSymInfo(xtal::BasicStructure const &prim);

  /// \brief Construct using factor group in given order
  PrimSymInfo(std::vector<xtal::SymOp> const &factor_group_elements,
              xtal::BasicStructure const &prim);

  /// \brief Construct using given factor group
  PrimSymInfo(std::shared_ptr<SymGroup const> const &_factor_group,
              BasicStructure const &prim);

  /// \brief Structure factor group
  std::shared_ptr<SymGroup const> factor_group;

  /// \brief Factor group elements without translations
  ///
  /// Notes:
  /// - Degenerate symmetry operations will not be added
  std::shared_ptr<SymGroup const> point_group;

  /// \brief Lattice point group
  std::shared_ptr<SymGroup const> lattice_point_group;

  /// \brief Describes how xtal::UnitCellCoord permute under symmetry
  ///
  /// Usage:
  /// - xtal::UnitCellCoord integral_site_coordinate_after =
  /// sym::copy_apply(unitcellcoord_symgroup_rep[factor_group_index],
  /// integral_site_coordinate_before)
  sym_info::UnitCellCoordSymGroupRep unitcellcoord_symgroup_rep;

  /// \brief True if any occupation DoF
  bool has_occupation_dofs;

  /// \brief True if any permutation in occ_symgroup_rep is non-trivial
  bool has_aniso_occs;

  /// \brief Permutations describe occupant index transformation under symmetry
  ///
  /// Usage:
  /// \code
  /// Index occupant_index_after = occ_symgroup_rep
  ///                                  .at(group_element_index)
  ///                                  .at(sublattice_index_before)
  ///                                  .at(occupant_index_before);
  /// \endcode
  ///
  /// Note:
  /// - This describes cases such as discrete molecular orientations or
  /// occupants with spin where a symmetry operation may transform one discrete
  /// occupant into another *before* permutating among sites.
  sym_info::OccSymGroupRep occ_symgroup_rep;

  /// \brief Permutations describe atom position index transformation under
  /// symmetry
  ///
  /// Usage:
  /// \code
  /// Index atom_position_index_after = atom_position_symgroup_rep
  ///                                       .at(group_element_index)
  ///                                       .at(sublattice_index_before)
  ///                                       .at(occupant_index_before)
  ///                                       .at(atom_position_index_before);
  /// \endcode
  ///
  /// Note:
  /// - This describes symmetry operations transforming molecules, resulting in
  /// permutation of atoms among the symmetrically equivalent atom positions in
  /// the molecule.
  sym_info::AtomPositionSymGroupRep atom_position_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// \code
  /// Eigen::MatrixXd const &M = local_dof_symgroup_rep
  ///                                .at(dof_type)
  ///                                .at(group_element_index)
  ///                                .at(sublattice_index_before);
  /// Eigen::MatrixXd sublattice_local_dof_values_after =
  ///     M * sublattice_local_dof_values_before;
  /// \endcode
  ///
  /// Note:
  /// - For each group element there is one matrix representation per sublattice
  /// - Local DoF values transform using these symrep matrices *before*
  ///   permuting among sites.
  /// - If has_occupation_dofs==true, a matrix rep for transforming occupation
  ///   indicator variables is included with key "occ".
  ///
  std::map<DoFKey, sym_info::LocalDoFSymGroupRep> local_dof_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// \code
  /// Eigen::MatrixXd const &M = global_dof_symgroup_rep
  ///                                .at(dof_type)
  ///                                .at(group_element_index);
  /// Eigen::VectorXd global_dof_values_after = M * global_dof_values_before;
  /// \endcode
  ///
  std::map<DoFKey, sym_info::GlobalDoFSymGroupRep> global_dof_symgroup_rep;
};

}  // namespace config
}  // namespace CASM

#endif
