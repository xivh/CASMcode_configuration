#ifndef CASM_occ_events_OccSystem
#define CASM_occ_events_OccSystem

#include <memory>
#include <string>
#include <vector>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class UnitCellCoord;
}  // namespace xtal

namespace occ_events {

struct OccPosition;

/// \brief Defines the system for OccPosition / OccTrajectory / OccEvent
struct OccSystem {
  OccSystem(std::vector<std::string> const &_resevoir_components,
            std::shared_ptr<xtal::BasicStructure const> const &_prim);

  /// \brief Names of the occupants in the resevoir
  ///
  /// Notes:
  /// - OccPosition::occupant_index is an index into this
  ///   when `is_in_resevoir==true`
  std::vector<std::string> resevoir_components;

  /// \brief The prim
  std::shared_ptr<xtal::BasicStructure const> prim;

  OccPosition make_resevoir_position(std::string molecule_name) const;

  OccPosition make_molecule_position(
      xtal::UnitCellCoord const &integral_site_coordinate,
      std::string molecule_name) const;

  OccPosition make_atomic_component_position(
      xtal::UnitCellCoord const &integral_site_coordinate,
      std::string molecule_name, Index atom_position_index) const;

  std::string get_name(OccPosition const &occ_position) const;

  Eigen::Vector3d get_cartesian_coordinate(
      OccPosition const &occ_position) const;
};

}  // namespace occ_events
}  // namespace CASM

#endif
