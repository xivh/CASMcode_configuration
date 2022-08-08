#include "casm/configuration/occ_events/OccTrajectory.hh"

#include "casm/configuration/occ_events/OccEventRep.hh"

namespace CASM {
namespace occ_events {

OccTrajectory::OccTrajectory() {}

OccTrajectory::OccTrajectory(std::initializer_list<OccPosition> _position)
    : position(_position) {}

OccTrajectory::OccTrajectory(std::vector<OccPosition> const &_position)
    : position(_position) {}

OccTrajectory::OccTrajectory(std::vector<OccPosition> &&_position)
    : position(std::move(_position)) {}

/// \brief Translate the OccEvent by a UnitCell translation
OccTrajectory &OccTrajectory::operator+=(xtal::UnitCell trans) {
  for (auto it = this->position.begin(); it != this->position.end(); ++it) {
    *it += trans;
  }
  return *this;
}

/// Compare OccTrajectory
bool OccTrajectory::operator<(OccTrajectory const &B) const {
  if (this->position.size() != B.position.size()) {
    return this->position.size() < B.position.size();
  }
  return this->position < B.position;
}

/// \brief Apply SymOp to OccTrajectory
OccTrajectory &apply(OccEventRep const &rep, OccTrajectory &occ_trajectory) {
  for (auto &occ_position : occ_trajectory.position) {
    apply(rep, occ_position);
  }
  return occ_trajectory;
}

/// \brief Apply SymOp to OccTrajectory
OccTrajectory copy_apply(OccEventRep const &rep, OccTrajectory occ_trajectory) {
  apply(rep, occ_trajectory);
  return occ_trajectory;
}

}  // namespace occ_events
}  // namespace CASM
