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

}  // namespace sym_info
}  // namespace CASM

#endif
