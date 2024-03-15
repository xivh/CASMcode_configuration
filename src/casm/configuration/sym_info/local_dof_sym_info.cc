#include "casm/configuration/sym_info/local_dof_sym_info.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace sym_info {

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
///
std::map<DoFKey, LocalDoFSymGroupRep> make_local_dof_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure) {
  auto const &basis = basicstructure.basis();
  double xtal_tol = basicstructure.lattice().tol();

  // sitemap[i]->std::vector<UnitCellCoord>:
  // - the vector of sites that i-th op maps prim basis sites onto
  // - used for constructing the site-related symreps
  std::vector<std::vector<UnitCellCoord>> sitemap;
  for (SymOp const &op : group_elements) {
    sitemap.push_back(xtal::symop_site_map(op, basicstructure));
  }

  std::map<DoFKey, LocalDoFSymGroupRep> local_dof_symgroup_rep;
  for (std::string dof_key : xtal::continuous_local_dof_types(basicstructure)) {
    LocalDoFSymGroupRep group_rep;
    Index op_index = 0;
    for (auto const &op : group_elements) {
      LocalDoFSymOpRep op_rep(basis.size());
      for (Index from_b = 0; from_b < basis.size(); ++from_b) {
        if (!basis[from_b].has_dof(dof_key)) continue;

        xtal::SiteDoFSet const &dof_from = basis[from_b].dof(dof_key);
        Index to_b = sitemap[op_index][from_b].sublattice();
        xtal::SiteDoFSet const &dof_to = basis[to_b].dof(dof_key);

        // want to check if:
        //   copy_apply(op, dof_from.basis()) = dof_to.basis() * U
        //               transformed_dof_from = dof_to.basis() * U

        xtal::SiteDoFSet transformed_dof_from = sym::copy_apply(op, dof_from);
        xtal::SiteDoFSetIsEquivalent_f dof_equals(dof_to, xtal_tol);

        if (!dof_equals(transformed_dof_from)) {
          std::stringstream msg;
          msg << "Error in make_local_dof_symgroup_rep: Local DoF \"" << dof_key
              << "\" that originally identified as equivalent cannot be mapped "
                 "by symmetry.";
          throw std::runtime_error(msg.str());
        }

        // B_to * x_to = op * B_from * x_from
        // x_to = M * x_from
        // ->
        // B_to * M * x_from = op * B_from * x_from
        // M = B_to.solve(op * B_from)

        Eigen::MatrixXd M;
        try {
          M = dof_to.basis().colPivHouseholderQr().solve(
              transformed_dof_from.basis());
          if (!(M.transpose() * M).eval().isIdentity(xtal_tol)) {
            throw std::runtime_error(
                "Cannot find orthogonal symmetry representation!");
          }
        } catch (std::runtime_error &e) {
          std::stringstream msg;
          msg << "Error in make_local_dof_symgroup_rep: Local DoF \"" << dof_key
              << "\" basis change representation construction failed.";
          throw std::runtime_error(msg.str());
        }
        op_rep[from_b] = M;
      }
      group_rep.push_back(op_rep);
      ++op_index;
    }
    local_dof_symgroup_rep.emplace(dof_key, group_rep);
  }
  return local_dof_symgroup_rep;
}

}  // namespace sym_info
}  // namespace CASM
