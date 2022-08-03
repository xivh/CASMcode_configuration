#ifndef CASM_sym_info_unitcellcoord
#define CASM_sym_info_unitcellcoord

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {
namespace sym_info {

UnitCellCoordSymGroupRep make_unitcellcoord_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure);

}
}  // namespace CASM

#endif
