#ifndef CASM_config_PrimSymInfo
#define CASM_config_PrimSymInfo

#include <map>

#include "casm/configuration/definitions.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace config {

/// \brief Data structure describing how prim DoF values transform under
/// application of symmetry
struct PrimSymInfo {
  PrimSymInfo(BasicStructure const &basicstructure);

  /// \brief Structure factor group
  std::shared_ptr<SymGroup const> factor_group;

  /// \brief Factor group elements without translations
  ///
  /// Notes:
  /// - Degenerate symmetry operations will not be added
  std::shared_ptr<SymGroup const> point_group;

  /// \brief Describes how xtal::UnitCellCoord permute under symmetry
  ///
  /// Usage:
  /// - xtal::UnitCellCoord integral_site_coordinate_after =
  /// sym::copy_apply(basis_permutation_symgroup_rep[factor_group_index],
  /// integral_site_coordinate_before)
  BasisPermutationSymGroupRep basis_permutation_symgroup_rep;

  /// \brief True if any occupation DoF
  bool has_occupation_dofs;

  /// \brief True if any permutation in occ_symgroup_rep is non-trivial
  bool has_aniso_occs;

  /// \brief Permutations describe occupant index transformation under symmetry
  ///
  /// Usage:
  /// - Index occupant_index_after =
  /// occ_symgroup_rep[factor_group_index][sublattice_index_before][occupant_index_before]
  /// - This rep is used by:
  ///   - `Configuration copy_apply(SupercellSymOp const &op, Configuration
  ///   configuration)`
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
  /// - Index atom_position_index_after =
  /// atom_position_symgroup_rep[factor_group_index][sublattice_index_before][occupant_index_before]
  ///
  /// Note:
  /// - This describes symmetry operations transforming molecules, resulting in
  /// permutation of atoms among the symmetrically equivalent atom positions in
  /// the molecule.
  AtomPositionSymGroupRep atom_position_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// -
  /// local_dof_symgroup_rep[dof_type][factor_group_index][sublattice_index_before]
  ///     -> Eigen::MatrixXd
  /// - For each group element there is one matrix representation per sublattice
  /// - This rep is used by:
  ///   - `Configuration copy_apply(SupercellSymOp const &op, Configuration
  ///   configuration)`
  ///
  /// Note:
  /// - Local DoF values transform using these symrep matrices *before*
  /// permuting among sites
  std::map<DoFKey, LocalDoFSymGroupRep> local_dof_symgroup_rep;

  /// \brief Matrices describe local DoF value transformation under symmetry
  ///
  /// Usage:
  /// - global_dof_symgroup_rep[dof_type][factor_group_index]
  ///       -> Eigen::MatrixXd
  /// - This rep is used by:
  ///   - `Configuration copy_apply(SupercellSymOp const &op, Configuration
  ///   configuration)`
  std::map<DoFKey, GlobalDoFSymGroupRep> global_dof_symgroup_rep;
};

}  // namespace config
}  // namespace CASM

#endif
