#ifndef CASM_sym_info_local_dof
#define CASM_sym_info_local_dof

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace sym_info {

/// \brief Matrices describe local DoF value transformation under symmetry
std::map<DoFKey, LocalDoFSymGroupRep> make_local_dof_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure);

}  // namespace sym_info
}  // namespace CASM

#endif
