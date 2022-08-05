#include "casm/configuration/occ_events/OccSystem.hh"

#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace occ_events {

OccSystem::OccSystem(std::vector<std::string> const &_resevoir_components,
                     std::shared_ptr<xtal::BasicStructure const> const &_prim)
    : resevoir_components(_resevoir_components), prim(_prim) {}

/// Make an OccPosition that indicates occupant in the resevoir
OccPosition OccSystem::make_resevoir_position(std::string molecule_name) const {
  Index occupant_index = find_index(this->resevoir_components, molecule_name);
  if (occupant_index < 0 ||
      occupant_index >= this->resevoir_components.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_resevoir_position: Invalid OccPosition "
        "molecule_name");
  }
  return OccPosition{true, false, xtal::UnitCellCoord{}, occupant_index, -1};
}

/// Make an OccPosition for a molecule (single atom or multi-atom)
OccPosition OccSystem::make_molecule_position(
    xtal::UnitCellCoord const &integral_site_coordinate,
    std::string molecule_name) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "integral_site_coordinate");
  }

  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  Index occupant_index = find_index_if(
      occupant_dof,
      [&](xtal::Molecule const &mol) { return mol.name() == molecule_name; });
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid molecule_name");
  }

  return OccPosition{false, false, integral_site_coordinate, occupant_index,
                     -1};
}

/// Make an OccPosition for a single atom in a multi-atom molecule
OccPosition OccSystem::make_atomic_component_position(
    xtal::UnitCellCoord const &integral_site_coordinate,
    std::string molecule_name, Index atom_position_index) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_component_position: Invalid "
        "integral_site_coordinate");
  }

  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  Index occupant_index = find_index_if(
      occupant_dof,
      [&](xtal::Molecule const &mol) { return mol.name() == molecule_name; });
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_component_position: Invalid "
        "molecule_name");
  }

  xtal::Molecule const &mol = occupant_dof[occupant_index];
  if (atom_position_index < 0 || atom_position_index >= mol.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_component_position: Invalid "
        "atom_position_index");
  }

  return OccPosition{false, true, integral_site_coordinate, occupant_index,
                     atom_position_index};
}

std::string OccSystem::get_name(OccPosition const &occ_position) const {
  auto const &pos = occ_position;
  if (pos.is_in_resevoir) {
    if (pos.occupant_index < 0 ||
        pos.occupant_index >= this->resevoir_components.size()) {
      throw std::runtime_error(
          "Error in OccSystem::get_name: Invalid OccPosition occupant_index in "
          "resevoir");
    }
    return this->resevoir_components[pos.occupant_index];
  }
  Index b = pos.integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_name: Invalid OccPosition sublattice");
  }
  xtal::Site const &site = this->prim->basis()[b];
  if (pos.occupant_index < 0 ||
      pos.occupant_index >= site.occupant_dof().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_name: Invalid OccPosition occupant_index");
  }
  xtal::Molecule const &mol = site.occupant_dof()[pos.occupant_index];
  if (!pos.is_atom) {
    return mol.name();
  }
  if (pos.atom_position_index < 0 || pos.atom_position_index >= mol.size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_name: Invalid OccPosition "
        "atom_position_index");
  }
  return mol.atom(pos.atom_position_index).name();
}

Eigen::Vector3d OccSystem::get_cartesian_coordinate(
    OccPosition const &occ_position) const {
  auto const &pos = occ_position;
  if (pos.is_in_resevoir) {
    throw std::runtime_error(
        "Error in OccSystem::get_cartesian_coordinate: OccPosition is in "
        "resevoir");
  }

  Index b = pos.integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_cartesian_coordinate: Invalid OccPosition "
        "sublattice");
  }
  xtal::Site const &site = this->prim->basis()[b];
  if (pos.occupant_index < 0 ||
      pos.occupant_index >= site.occupant_dof().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_cartesian_coordinate: Invalid OccPosition "
        "occupant_index");
  }
  Eigen::Vector3d coord = site.const_cart();
  if (!pos.is_atom) {
    return coord;
  }
  xtal::Molecule const &mol = site.occupant_dof()[pos.occupant_index];
  if (pos.atom_position_index < 0 || pos.atom_position_index >= mol.size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_cartesian_coordinate: Invalid OccPosition "
        "atom_position_index");
  }
  coord += mol.atom(pos.atom_position_index).cart();
  return coord;
}

}  // namespace occ_events
}  // namespace CASM
