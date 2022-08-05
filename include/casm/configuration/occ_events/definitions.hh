// The `casm/configuration/occ_events` module supports
// construction of orbits of KMC-type occupation events
//
// Primarily, this purpose of this module is to provide:
// - OccEvent, a vector of OccTrajectory
// - OccTrajectory, a vector of OccPosition
// - OccPosition, specifies a molecule on a site, an atom in a
//   molecule on a site, or a molecule in a resovoir
// - Methods for finding orbits of OccEvent and the symmetry
//   of OccEvent
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_crystallography
// - CASMcode_configuration/group
// - CASMcode_configuration/sym_info
// - CASMcode_configuration/clusterography

#ifndef CASM_occ_events_definitions
#define CASM_occ_events_definitions

#include "casm/configuration/clusterography/definitions.hh"
#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class Site;
struct SymOp;
class UnitCell;
class UnitCellCoord;
struct UnitCellCoordRep;
}  // namespace xtal

namespace group {
template <typename ElementType>
struct Group;
}

namespace occ_events {

typedef long Index;

/// \brief A group::Group of xtal::SymOp
typedef group::Group<xtal::SymOp> SymGroup;

class OccEvent;
class OccEventInvariants;
struct OccEventRep;
struct OccPosition;
struct OccSystem;
struct OccTrajectory;

}  // namespace occ_events
}  // namespace CASM

#endif
