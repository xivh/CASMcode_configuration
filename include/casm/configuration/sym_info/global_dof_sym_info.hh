#ifndef CASM_sym_info_global_dof
#define CASM_sym_info_global_dof

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace sym_info {

/// \brief Matrices describe global DoF value transformation under symmetry
std::map<DoFKey, GlobalDoFSymGroupRep> make_global_dof_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure);

}  // namespace sym_info
}  // namespace CASM

#endif
