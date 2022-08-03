#include "casm/configuration/sym_info/factor_group.hh"

#include "casm/configuration/group/Group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SymTypeComparator.hh"

namespace CASM {
namespace sym_info {

/// \brief Generate prim factor group
///
/// Notes:
/// - Degenerate symmetry operations will not be added
std::shared_ptr<SymGroup const> make_factor_group(
    xtal::BasicStructure const &prim) {
  std::vector<SymOp> factor_group_elements = xtal::make_factor_group(prim);
  xtal::Lattice const &lattice = prim.lattice();
  double xtal_tol = lattice.tol();
  std::multiplies<SymOp> multiply_f;
  xtal::SymOpPeriodicCompare_f equal_to_f(lattice, xtal_tol);
  return std::make_shared<SymGroup>(
      group::make_group(factor_group_elements, multiply_f, equal_to_f));
}

/// \brief Generate prim point group
///
/// Notes:
/// - Degenerate symmetry operations will not be added
std::shared_ptr<SymGroup const> make_point_group(
    xtal::BasicStructure const &prim,
    std::shared_ptr<SymGroup const> const &factor_group) {
  xtal::Lattice const &lattice = prim.lattice();
  double xtal_tol = lattice.tol();
  std::multiplies<SymOp> multiply_f;
  xtal::SymOpPeriodicCompare_f equal_to_f(lattice, xtal_tol);

  std::vector<SymOp> pg_element;
  for (SymOp const &op : factor_group->element) {
    SymOp point_op(op.matrix, xtal::SymOpTranslationType::Zero(),
                   op.is_time_reversal_active);
    auto unary_f = [&](SymOp const &lhs) { return equal_to_f(lhs, point_op); };
    auto it = std::find_if(pg_element.begin(), pg_element.end(), unary_f);
    if (it == pg_element.end()) {
      pg_element.push_back(point_op);
    }
  }

  return std::make_shared<SymGroup>(
      group::make_group(pg_element, multiply_f, equal_to_f));
}

}  // namespace sym_info
}  // namespace CASM
