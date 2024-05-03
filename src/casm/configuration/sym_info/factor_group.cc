#include "casm/configuration/sym_info/factor_group.hh"

#include "casm/configuration/group/Group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/SymTypeComparator.hh"

namespace CASM {
namespace sym_info {

ConjugacyClassCompare::ConjugacyClassCompare(double _tol) : tol(_tol) {}

bool ConjugacyClassCompare::operator()(const map_type &A,
                                       const map_type &B) const {
  return float_lexicographical_compare(A.begin()->first, B.begin()->first, tol);
}

/// \brief Construct a SymGroup from the group elements, sorting by class
/// and operation
std::shared_ptr<SymGroup const> make_symgroup(
    std::vector<SymOp> const &_elements, xtal::Lattice const &lattice) {
  std::vector<SymOp> elements{_elements};

  // this is called `sort_factor_group`, but can sort any xtal::SymOp
  sort_factor_group(elements, lattice);

  double xtal_tol = lattice.tol();
  std::multiplies<SymOp> multiply_f;
  xtal::SymOpPeriodicCompare_f equal_to_f(lattice, xtal_tol);

  // make a `tmp` group, so we can use the multiplication table
  SymGroup tmp = group::make_group(elements, multiply_f, equal_to_f);

  // use the group with multiplication table to make conjugacy classes
  std::vector<std::vector<Index>> conjugacy_classes =
      make_conjugacy_classes(tmp);

  // insert elements into each conjugacy class to sort
  // elements in the class (using a std::map),
  // then sort the classes by the first element in the
  // class (using a std::set)
  xtal::SymOpSortKeyCompare op_compare(xtal_tol);
  ConjugacyClassCompare cclass_compare(xtal_tol);
  std::set<ConjugacyClassCompare::map_type, ConjugacyClassCompare> sorter(
      cclass_compare);
  for (int i = 0; i < conjugacy_classes.size(); ++i) {
    ConjugacyClassCompare::map_type cclass(op_compare);
    for (int j = 0; j < conjugacy_classes[i].size(); ++j) {
      const SymOp &op = tmp.element.at(conjugacy_classes[i][j]);
      cclass.emplace(make_symop_sort_key(op, lattice), op);
    }
    sorter.emplace(std::move(cclass));
  }

  // copy the elements back out now that they are sorted by class
  std::vector<xtal::SymOp> sorted_elements;
  for (auto const &cclass : sorter) {
    for (auto const &pair : cclass) {
      sorted_elements.push_back(pair.second);
    }
  }

  // now make the group with elements sorted by class
  return std::make_shared<SymGroup>(
      group::make_group(sorted_elements, multiply_f, equal_to_f));
}

/// \brief Construct a SymGroup from the group elements,
///     without sorting by class or operation
std::shared_ptr<SymGroup const> make_symgroup_without_sorting(
    std::vector<SymOp> const &elements, xtal::Lattice const &lattice) {
  double xtal_tol = lattice.tol();
  std::multiplies<SymOp> multiply_f;
  xtal::SymOpPeriodicCompare_f equal_to_f(lattice, xtal_tol);
  return std::make_shared<SymGroup>(
      group::make_group(elements, multiply_f, equal_to_f));
}

/// \brief Generate prim factor group
///
/// Notes:
/// - Result is sorted by class, with classes sorted by symop
/// - Uses lattice tol for comparison
std::shared_ptr<SymGroup const> make_factor_group(
    xtal::BasicStructure const &prim) {
  std::vector<SymOp> elements = xtal::make_factor_group(prim);
  return make_symgroup(elements, prim.lattice());
}

/// \brief Use prim factor group
std::shared_ptr<SymGroup const> use_factor_group(
    std::vector<xtal::SymOp> const &factor_group_elements,
    xtal::BasicStructure const &prim) {
  return make_symgroup_without_sorting(factor_group_elements, prim.lattice());
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

/// \brief Generate lattice point group
std::shared_ptr<SymGroup const> make_lattice_point_group(
    xtal::Lattice const &lattice) {
  // these SymOp are not sorted
  std::vector<SymOp> elements = xtal::make_point_group(lattice);
  return make_symgroup(elements, lattice);
}

}  // namespace sym_info
}  // namespace CASM
