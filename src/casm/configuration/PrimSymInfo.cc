#include "casm/configuration/PrimSymInfo.hh"

#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/configuration/sym_info/global_dof_sym_info.hh"
#include "casm/configuration/sym_info/local_dof_sym_info.hh"
#include "casm/configuration/sym_info/occ_sym_info.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace config {

/// \brief Constructor
///
/// Notes:
/// - Uses `sym_info::make_factor_group(prim);`
PrimSymInfo::PrimSymInfo(BasicStructure const &prim)
    : PrimSymInfo(sym_info::make_factor_group(prim), prim) {}

/// \brief Construct using factor group in given order
///
/// Notes:
/// - Uses `sym_info::use_factor_group(factor_group_elements, prim);`
PrimSymInfo::PrimSymInfo(std::vector<xtal::SymOp> const &factor_group_elements,
                         BasicStructure const &prim)
    : PrimSymInfo(sym_info::use_factor_group(factor_group_elements, prim),
                  prim) {}

/// \brief Construct using given factor group
PrimSymInfo::PrimSymInfo(std::shared_ptr<SymGroup const> const &_factor_group,
                         BasicStructure const &prim)
    : factor_group(_factor_group) {
  using namespace sym_info;

  this->point_group = make_point_group(prim, this->factor_group);
  this->lattice_point_group = make_lattice_point_group(prim.lattice());

  this->unitcellcoord_symgroup_rep =
      make_unitcellcoord_symgroup_rep(this->factor_group->element, prim);

  OccSymInfo occ_sym_info(this->factor_group->element, prim);
  this->has_occupation_dofs = occ_sym_info.has_occupation_dofs;
  this->has_aniso_occs = occ_sym_info.has_aniso_occs;
  this->occ_symgroup_rep = occ_sym_info.occ_symgroup_rep;
  this->atom_position_symgroup_rep = occ_sym_info.atom_position_symgroup_rep;

  this->local_dof_symgroup_rep =
      make_local_dof_symgroup_rep(this->factor_group->element, prim);

  if (this->has_occupation_dofs) {
    local_dof_symgroup_rep["occ"] = sym_info::LocalDoFSymGroupRep();
    for (sym_info::OccSymOpRep const &occ_symop_rep : this->occ_symgroup_rep) {
      std::vector<Eigen::MatrixXd> occ_matrix_rep;
      for (sym_info::Permutation const &perm : occ_symop_rep) {
        occ_matrix_rep.push_back(sym_info::as_matrix(sym_info::inverse(perm)));
      }
      local_dof_symgroup_rep.at("occ").push_back(occ_matrix_rep);
    }
  }

  this->global_dof_symgroup_rep =
      make_global_dof_symgroup_rep(this->factor_group->element, prim);
}

}  // namespace config
}  // namespace CASM
