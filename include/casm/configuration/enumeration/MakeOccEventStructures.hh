#ifndef CASM_config_enum_MakeOccEventStructures
#define CASM_config_enum_MakeOccEventStructures

#include <memory>
#include <string>
#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
}

namespace occ_events {
class OccEvent;
struct OccSystem;
}  // namespace occ_events

namespace config {

struct Configuration;

/// \brief Generate xtal::SimpleStructure, properly aligned for NEB
///     calculations
///
/// Notes:
/// - Only occupation and strain DoF are allowed
/// - Vacancies are not included
/// - Throws if Multi-atom molecules
/// - Throws if OccTrajectory are not size 2 (i.e. intermediate event states)
/// - Atom types are sorted according to
/// occ_events::OccSystem::chemical_name_list, with
///   atoms involved in the event at the beginning of each section
class MakeOccEventStructures {
 public:
  /// \brief MakeOccEventStructures constructor
  MakeOccEventStructures(
      Configuration const &configuration, occ_events::OccEvent const &occ_event,
      std::shared_ptr<occ_events::OccSystem const> const &system,
      bool skip_event_occupants = false);

  /// \brief Construct interpolated structure
  xtal::SimpleStructure operator()(double interpolation_factor) const;

 private:
  /// Ideal lattice
  Eigen::Matrix3d m_ideal_lat_column_mat;

  /// Atom names
  std::vector<std::string> m_atom_names;

  /// Cartesian coordinates at initial structure
  std::vector<Eigen::Vector3d> m_coords_init;

  /// Displacement from initial to final structure
  /// (atoms not involved in event will have {0., 0., 0.})
  std::vector<Eigen::Vector3d> m_event_disp;

  /// Configuration deformation gradient
  Eigen::Matrix3d m_F;
};

}  // namespace config
}  // namespace CASM

#endif
