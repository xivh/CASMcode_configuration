#ifndef CASM_occ_events_OccEventRep
#define CASM_occ_events_OccEventRep

#include "casm/configuration/occ_events/definitions.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace occ_events {

/// \brief Transform OccPosition / OccTrajectory / OccEvent
struct OccEventRep {
  OccEventRep(xtal::UnitCellCoordRep const &_unitcellcoord_rep,
              sym_info::OccSymOpRep const &_occupant_rep,
              sym_info::AtomPositionSymOpRep const &_atom_position_rep)
      : unitcellcoord_rep(_unitcellcoord_rep),
        occupant_rep(_occupant_rep),
        atom_position_rep(_atom_position_rep) {}

  xtal::UnitCellCoordRep unitcellcoord_rep;
  sym_info::OccSymOpRep occupant_rep;
  sym_info::AtomPositionSymOpRep atom_position_rep;
};

std::vector<OccEventRep> make_occevent_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &prim);

std::vector<OccEventRep> make_occevent_symgroup_rep(
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
    std::vector<sym_info::OccSymOpRep> const &occ_symgroup_rep,
    std::vector<sym_info::AtomPositionSymOpRep> const
        &atom_position_symgroup_rep);

}  // namespace occ_events
}  // namespace CASM

#endif
