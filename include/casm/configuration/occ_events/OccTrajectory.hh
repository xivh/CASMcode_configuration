#ifndef CASM_occ_events_OccTrajectory
#define CASM_occ_events_OccTrajectory

#include <vector>

#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/sym_info/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace occ_events {

/// A single occupant trajectory as a vector of OccPosition
struct OccTrajectory
    : public Comparisons<xtal::Translatable<CRTPBase<OccTrajectory>>> {
  /// A vector of OccPosition
  std::vector<OccPosition> position;

  /// Translate the OccTrajectory by a UnitCell translation
  OccTrajectory &operator+=(xtal::UnitCell trans);

  /// Compare OccTrajectory
  bool operator<(OccTrajectory const &B) const;
};

/// \brief Apply SymOp to OccTrajectory
OccTrajectory &apply(OccEventRep const &rep, OccTrajectory &occ_trajectory);

/// \brief Apply SymOp to OccTrajectory
OccTrajectory copy_apply(OccEventRep const &rep, OccTrajectory occ_trajectory);

}  // namespace occ_events
}  // namespace CASM

#endif
