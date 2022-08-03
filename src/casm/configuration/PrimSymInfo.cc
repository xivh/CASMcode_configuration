#include "casm/configuration/PrimSymInfo.hh"

#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/configuration/sym_info/global_dof_sym_info.hh"
#include "casm/configuration/sym_info/local_dof_sym_info.hh"
#include "casm/configuration/sym_info/occ_sym_info.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"

namespace CASM {
namespace config {

/// \brief Constructor
///
/// Notes:
/// - `std::vector<SymOp> fg_element = xtal::make_factor_group(basicstructure);`
/// - Uses basicstructure.lattice().tol() for symmetry finding
PrimSymInfo::PrimSymInfo(BasicStructure const &prim) {
  using namespace sym_info;

  this->factor_group = make_factor_group(prim);
  this->point_group = make_point_group(prim, this->factor_group);

  this->unitcellcoord_symgroup_rep =
      make_unitcellcoord_symgroup_rep(this->factor_group->element, prim);

  OccSymInfo occ_sym_info(this->factor_group->element, prim);
  this->has_occupation_dofs = occ_sym_info.has_occupation_dofs;
  this->has_aniso_occs = occ_sym_info.has_aniso_occs;
  this->occ_symgroup_rep = occ_sym_info.occ_symgroup_rep;
  this->atom_position_symgroup_rep = occ_sym_info.atom_position_symgroup_rep;

  this->local_dof_symgroup_rep =
      make_local_dof_symgroup_rep(this->factor_group->element, prim);
  this->global_dof_symgroup_rep =
      make_global_dof_symgroup_rep(this->factor_group->element, prim);
}

}  // namespace config
}  // namespace CASM
