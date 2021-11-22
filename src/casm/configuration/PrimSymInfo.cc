#include <algorithm>
#include "casm/configuration/PrimSymInfo.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SymTypeComparator.hh"

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

  // Construct std::shared_ptr<Group const> factor_group;
  std::vector<SymOp> fg_element = xtal::make_factor_group(basicstructure);
  factor_group = std::make_shared<SymGroup>(
      make_group(fg_element, multiply_f, equal_to_f));

  // Construct std::shared_ptr<Group const> point_group;
  // Notes:
  // - Degenerate symmetry operations will not be added
  std::vector<SymOp> pg_element;
  for (SymOp const &op : fg_element) {
    SymOp point_op(op.matrix, xtal::SymOpTranslationType::Zero(),
                   op.is_time_reversal_active);
    auto unary_f = [&](SymOp const &lhs) {
      return equal_to_f(lhs, point_op);
    };
    auto it = std::find_if(pg_element.begin(), pg_element.end(), unary_f);
    if (it == pg_element.end()) {
      pg_element.push_back(point_op);
    }
  }
  point_group = std::make_shared<SymGroup>(
      make_group(pg_element, multiply_f, equal_to_f));

  // BasisPermutationSymGroupRep basis_permutation_symgroup_rep;

  // bool has_aniso_occs;

  // OccSymGroupRep occ_symgroup_rep;

  // std::map<DoFKey, LocalDoFSymGroupRep>  local_dof_symgroup_rep;

  // std::map<DoFKey, GlobalDoFSymGroupRep> global_dof_symgroup_rep;
}

}  // namespace config
}  // namespace CASM
