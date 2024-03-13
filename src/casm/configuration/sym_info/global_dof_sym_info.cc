#include "casm/configuration/sym_info/global_dof_sym_info.hh"

#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace sym_info {

/// \brief Matrices describe global DoF value transformation under symmetry
///
/// Usage:
/// \code
/// Eigen::MatrixXd const &M = global_dof_symgroup_rep
///                                .at(dof_type)
///                                .at(group_element_index);
/// Eigen::VectorXd global_dof_values_after = M * global_dof_values_before;
/// \endcode
///
std::map<DoFKey, GlobalDoFSymGroupRep> make_global_dof_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure) {
  double xtal_tol = basicstructure.lattice().tol();

  std::map<DoFKey, GlobalDoFSymGroupRep> global_dof_symgroup_rep;
  for (auto const &pair : basicstructure.global_dofs()) {
    std::string const &dof_key = pair.first;
    xtal::DoFSet const &dof = pair.second;

    GlobalDoFSymGroupRep group_rep;
    for (auto const &op : group_elements) {
      xtal::DoFSetIsEquivalent_f dof_equals(dof, xtal_tol);
      xtal::DoFSet transformed_dof = sym::copy_apply(op, dof);
      if (!dof_equals(transformed_dof)) {
        std::stringstream msg;
        msg << "Error in CASM::config::Prim constructor: Global DoF \""
            << dof_key
            << "\" that originally identified as equivalent cannot be mapped "
               "by symmetry.";
        throw std::runtime_error(msg.str());
      }

      // B * x_after = op * B * x_before
      // x_after = M * x_before
      // ->
      // B * M * x_before = op * B * x_before
      // M = B.solve(op * B)

      Eigen::MatrixXd M;
      try {
        M = dof.basis().colPivHouseholderQr().solve(transformed_dof.basis());
        if (!(M.transpose() * M).eval().isIdentity(xtal_tol)) {
          throw std::runtime_error(
              "Cannot find orthogonal symmetry representation!");
        }
      } catch (std::runtime_error &e) {
        std::stringstream msg;
        msg << "Error in CASM::config::Prim constructor: Global DoF \""
            << dof_key << "\" basis change representation construction failed.";
        throw std::runtime_error(msg.str());
      }
      group_rep.push_back(M);
    }

    global_dof_symgroup_rep.emplace(dof_key, group_rep);
  }
  return global_dof_symgroup_rep;
}

}  // namespace sym_info
}  // namespace CASM
