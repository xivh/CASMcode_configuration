#include "casm/configuration/PrimSymInfo.hh"

#include <algorithm>

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/OccupantDoFIsEquivalent.hh"
#include "casm/crystallography/SymTypeComparator.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace config {

/// \brief Constructor
///
/// Notes:
/// - Uses basicstructure.lattice().tol() for symmetry finding
PrimSymInfo::PrimSymInfo(BasicStructure const &basicstructure) {
  // TODO: set the following

  xtal::Lattice const &lattice = basicstructure.lattice();
  double xtal_tol = lattice.tol();
  std::multiplies<SymOp> multiply_f;
  xtal::SymOpPeriodicCompare_f equal_to_f(lattice, xtal_tol);
  using basic_symmetry::make_group;

  // Construct (std::shared_ptr<Group const>) factor_group
  std::vector<SymOp> fg_element = xtal::make_factor_group(basicstructure);
  factor_group = std::make_shared<SymGroup>(
      make_group(fg_element, multiply_f, equal_to_f));

  // Construct (std::shared_ptr<Group const>) point_group
  // Notes:
  // - Degenerate symmetry operations will not be added
  std::vector<SymOp> pg_element;
  for (SymOp const &op : fg_element) {
    SymOp point_op(op.matrix, xtal::SymOpTranslationType::Zero(),
                   op.is_time_reversal_active);
    auto unary_f = [&](SymOp const &lhs) { return equal_to_f(lhs, point_op); };
    auto it = std::find_if(pg_element.begin(), pg_element.end(), unary_f);
    if (it == pg_element.end()) {
      pg_element.push_back(point_op);
    }
  }
  point_group = std::make_shared<SymGroup>(
      make_group(pg_element, multiply_f, equal_to_f));

  // sitemap[i]->std::vector<UnitCellCoord>:
  // - the vector of sites that i-th op maps prim basis sites onto
  // - used for constructing the site-related symreps
  Index op_index;
  std::vector<std::vector<UnitCellCoord>> sitemap;
  for (SymOp const &op : fg_element) {
    sitemap.push_back(xtal::symop_site_map(op, basicstructure));
  }

  // Construct (BasisPermutationSymGroupRep) basis_permutation_symgroup_rep
  op_index = 0;
  for (SymOp const &op : fg_element) {
    UnitCellCoordRep rep;
    rep.point_matrix = lround(cart2frac(op.matrix, lattice));
    for (UnitCellCoord const &site : sitemap[op_index]) {
      rep.sublattice_index.push_back(site.sublattice());
      rep.unitcell_indices.push_back(site.unitcell());
    }
    basis_permutation_symgroup_rep.push_back(rep);
    ++op_index;
  }

  // Construct (OccSymGroupRep) occ_symgroup_rep & set (bool) has_aniso_occs;
  auto const &basis = basicstructure.basis();
  op_index = 0;
  has_aniso_occs = false;
  for (SymOp const &op : fg_element) {
    // the rep is a vector of permutation (one per prim site), that are used
    // to generate the new occupant indices
    OccSymOpRep rep;
    for (Index b = 0; b < basis.size(); ++b) {
      // copy_aply(symop,dofref_from) = P.permute(dofref_to);
      auto const &dofref_to =
          basis[sitemap[op_index][b].sublattice()].occupant_dof();
      auto const &dofref_from = basis[b].occupant_dof();

      xtal::OccupantDoFIsEquivalent eq(dofref_from);

      if (eq(op, dofref_to)) {
        if (!eq.perm().is_identity()) {
          has_aniso_occs = true;
        }
        rep.push_back(inverse(eq.perm().perm_array()));
      } else {
        throw std::runtime_error(
            "In CASM::config::Prim constructor: Sites originally "
            "identified as equivalent cannot be mapped by symmetry.");
      }
    }
    occ_symgroup_rep.push_back(rep);
    ++op_index;
  }

  // std::map<DoFKey, LocalDoFSymGroupRep>  local_dof_symgroup_rep;
  for (std::string dof_key : xtal::continuous_local_dof_types(basicstructure)) {
    LocalDoFSymGroupRep group_rep;
    op_index = 0;
    for (auto const &op : fg_element) {
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
          msg << "Error in CASM::config::Prim constructor: Local DoF \""
              << dof_key
              << "\" that originally identified as equivalent cannot be mapped "
                 "by symmetry.";
          throw std::runtime_error(msg.str());
        }
        Eigen::MatrixXd basis_change_representation;
        try {
          basis_change_representation = xtal::dofset_transformation_matrix(
              dof_to.basis(), transformed_dof_from.basis(), xtal_tol);
        } catch (std::runtime_error &e) {
          std::stringstream msg;
          msg << "Error in CASM::config::Prim constructor: Local DoF \""
              << dof_key
              << "\" basis change representation construction failed.";
          throw std::runtime_error(msg.str());
        }
        op_rep[from_b] = basis_change_representation;
      }
      group_rep.push_back(op_rep);
      ++op_index;
    }
    local_dof_symgroup_rep.emplace(dof_key, group_rep);
  }

  // Construct (std::map<DoFKey, GlobalDoFSymGroupRep>) global_dof_symgroup_rep
  for (auto const &pair : basicstructure.global_dofs()) {
    std::string const &dof_key = pair.first;
    xtal::DoFSet const &dof = pair.second;

    GlobalDoFSymGroupRep group_rep;
    for (auto const &op : fg_element) {
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
      Eigen::MatrixXd basis_change_representation;
      try {
        basis_change_representation = xtal::dofset_transformation_matrix(
            dof.basis(), transformed_dof.basis(), xtal_tol);
      } catch (std::runtime_error &e) {
        std::stringstream msg;
        msg << "Error in CASM::config::Prim constructor: Global DoF \""
            << dof_key << "\" basis change representation construction failed.";
        throw std::runtime_error(msg.str());
      }
      group_rep.push_back(basis_change_representation);
    }

    global_dof_symgroup_rep.emplace(dof_key, group_rep);
  }
}

}  // namespace config
}  // namespace CASM
