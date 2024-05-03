#ifndef CASM_sym_info_factor_group
#define CASM_sym_info_factor_group

#include "casm/configuration/sym_info/definitions.hh"
#include "casm/crystallography/BasicStructureTools.hh"

namespace CASM {
namespace sym_info {

struct ConjugacyClassCompare {
  typedef std::map<xtal::symop_sort_key_type, xtal::SymOp,
                   xtal::SymOpSortKeyCompare>
      map_type;

  double tol;
  ConjugacyClassCompare(double _tol);

  bool operator()(const map_type &A, const map_type &B) const;
};

/// \brief Construct a SymGroup from the group elements, sorting by class
/// and operation
std::shared_ptr<SymGroup const> make_symgroup(
    std::vector<SymOp> const &elements, xtal::Lattice const &lattice);

/// \brief Construct a SymGroup from the group elements,
///     without sorting by class
std::shared_ptr<SymGroup const> make_symgroup_without_sorting(
    std::vector<SymOp> const &elements, xtal::Lattice const &lattice);

/// \brief Generate prim factor group
std::shared_ptr<SymGroup const> make_factor_group(
    xtal::BasicStructure const &prim);

/// \brief Use prim factor group
std::shared_ptr<SymGroup const> use_factor_group(
    std::vector<xtal::SymOp> const &factor_group_elements,
    xtal::BasicStructure const &prim);

/// \brief Generate prim point group
std::shared_ptr<SymGroup const> make_point_group(
    xtal::BasicStructure const &prim,
    std::shared_ptr<SymGroup const> const &factor_group);

/// \brief Generate lattice point group
std::shared_ptr<SymGroup const> make_lattice_point_group(
    xtal::Lattice const &lattice);

}  // namespace sym_info
}  // namespace CASM

#endif
