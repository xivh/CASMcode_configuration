#include "casm/configuration/occ_events/OccEventRep.hh"

#include "casm/configuration/sym_info/occ_sym_info.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {
namespace occ_events {

std::vector<OccEventRep> make_occevent_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &prim) {
  using namespace sym_info;
  std::vector<xtal::UnitCellCoordRep> unitcellcoord_symgroup_rep =
      make_unitcellcoord_symgroup_rep(group_elements, prim);
  OccSymInfo occ_sym_info(group_elements, prim);

  return make_occevent_symgroup_rep(unitcellcoord_symgroup_rep,
                                    occ_sym_info.occ_symgroup_rep,
                                    occ_sym_info.atom_position_symgroup_rep);
}

std::vector<OccEventRep> make_occevent_symgroup_rep(
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
    std::vector<sym_info::OccSymOpRep> const &occ_symgroup_rep,
    std::vector<sym_info::AtomPositionSymOpRep> const
        &atom_position_symgroup_rep) {
  std::vector<OccEventRep> occevent_symgroup_rep;
  for (Index i = 0; i < unitcellcoord_symgroup_rep.size(); ++i) {
    occevent_symgroup_rep.emplace_back(unitcellcoord_symgroup_rep[i],
                                       occ_symgroup_rep[i],
                                       atom_position_symgroup_rep[i]);
  }
  return occevent_symgroup_rep;
}

}  // namespace occ_events
}  // namespace CASM
