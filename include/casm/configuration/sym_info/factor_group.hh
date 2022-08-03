#ifndef CASM_sym_info_factor_group
#define CASM_sym_info_factor_group

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace sym_info {

/// \brief Generate prim factor group
std::shared_ptr<SymGroup const> make_factor_group(
    xtal::BasicStructure const &prim);

/// \brief Generate prim point group
std::shared_ptr<SymGroup const> make_point_group(
    xtal::BasicStructure const &prim,
    std::shared_ptr<SymGroup const> const &factor_group);

}  // namespace sym_info
}  // namespace CASM

#endif
