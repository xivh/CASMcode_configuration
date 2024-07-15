#include "casm/configuration/occ_events/OccSystem.hh"

#include <optional>

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/algorithm.hh"

// debug
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"

namespace CASM {
namespace occ_events {

OccSystem::OccSystem(std::shared_ptr<xtal::BasicStructure const> const &_prim,
                     std::vector<std::string> const &_chemical_name_list,
                     std::set<std::string> const &_vacancy_name_list)
    : prim(_prim),
      chemical_name_list(_chemical_name_list),
      vacancy_name_list(_vacancy_name_list),
      atom_name_list(make_atom_name_list(*prim)),
      orientation_name_list(make_orientation_name_list(*prim)) {
  Index occupant_index;

  for (auto const &chemical_name : this->chemical_name_list) {
    this->is_vacancy_list.push_back(vacancy_name_list.count(chemical_name));
  }

  // --- build is_indivisible_chemical_list ---
  {
    std::vector<std::optional<bool>> tmp_list(this->chemical_name_list.size());
    for (Index b = 0; b < prim->basis().size(); ++b) {
      Index occupant_index = 0;
      for (auto const &mol : prim->basis()[b].occupant_dof()) {
        Index chemical_index = find_index(this->chemical_name_list, mol.name());
        std::optional<bool> &ref = tmp_list[chemical_index];
        if (ref.has_value()) {
          if (*ref != mol.is_indivisible()) {
            std::stringstream ss;
            ss << "Error: inconsistent indivisible molecule flag (basis site "
               << b << ", occupant " << occupant_index << ")";
            throw std::runtime_error(ss.str());
          }
        } else {
          ref = mol.is_indivisible();
        }
        ++occupant_index;
      }
    }
    for (Index chemical_index = 0; chemical_index < chemical_name_list.size();
         ++chemical_index) {
      auto const &ref = tmp_list[chemical_index];
      if (!ref.has_value()) {
        std::stringstream ss;
        ss << "Error: molecule not found in prim: "
           << chemical_name_list[chemical_index];
        throw std::runtime_error(ss.str());
      }
      is_indivisible_chemical_list.push_back(*ref);
    }
  }

  // --- build occupant_index lookup tables ---

  for (Index b = 0; b < prim->basis().size(); ++b) {
    // build occupant / chemical_name_list index lookup tables
    std::vector<int> occ_to_chem;
    occupant_index = 0;
    for (auto const &mol : prim->basis()[b].occupant_dof()) {
      Index chemical_index = find_index(this->chemical_name_list, mol.name());
      occ_to_chem.push_back(chemical_index);
      ++occupant_index;
    }
    occupant_to_chemical_index.push_back(occ_to_chem);

    // build occupant / orientation_name_list lookup tables
    std::vector<int> occ_to_orient;
    std::vector<int> orient_to_occ(this->orientation_name_list.size(), -1);
    occupant_index = 0;
    for (auto const &name : prim->unique_names()[b]) {
      Index orientation_index = find_index(this->orientation_name_list, name);

      occ_to_orient.push_back(orientation_index);
      orient_to_occ[orientation_index] = occupant_index;
      ++occupant_index;
    }
    occupant_to_orientation_index.push_back(occ_to_orient);
    orientation_to_occupant_index.push_back(orient_to_occ);
  }

  // --- build atom_position lookup table ---

  atom_position_to_name_index.resize(prim->basis().size());
  Index b = 0;
  for (auto const &site : prim->basis()) {
    atom_position_to_name_index[b].resize(site.occupant_dof().size());
    Index occupant_index = 0;
    for (auto const &mol : site.occupant_dof()) {
      for (auto const &atom : mol.atoms()) {
        Index atom_name_index = find_index(this->atom_name_list, atom.name());
        atom_position_to_name_index[b][occupant_index].push_back(
            atom_name_index);
      }
      ++occupant_index;
    }
    ++b;
  }
}

OccPosition OccSystem::make_molecule_position(
    xtal::UnitCellCoord const &integral_site_coordinate,
    Index occupant_index) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "integral_site_coordinate");
  }
  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid occupant_index");
  }

  return OccPosition{false, false, integral_site_coordinate, occupant_index,
                     -1};
}

OccPosition OccSystem::make_atom_position(
    xtal::UnitCellCoord const &integral_site_coordinate, Index occupant_index,
    Index atom_position_index) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "integral_site_coordinate");
  }
  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid occupant_index");
  }
  xtal::Molecule const &mol = occupant_dof[occupant_index];
  if (atom_position_index < 0 || atom_position_index >= mol.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "atom_position_index");
  }

  return OccPosition{false, true, integral_site_coordinate, occupant_index,
                     atom_position_index};
}

/// Make an OccPosition for a molecule (single atom or multi-atom)
OccPosition OccSystem::make_molecule_position(
    xtal::UnitCellCoord const &integral_site_coordinate,
    std::string orientation_name) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid "
        "integral_site_coordinate");
  }

  if (this->prim->unique_names().size() != this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: orientation_name & basis "
        "size mismatch");
  }
  if (this->prim->unique_names()[b].size() !=
      this->prim->basis()[b].occupant_dof().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: orientation_name & basis "
        "occupant_dof size mismatch");
  }

  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  Index occupant_index =
      find_index(this->prim->unique_names()[b], orientation_name);
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_position: Invalid orientation_name");
  }

  return OccPosition{false, false, integral_site_coordinate, occupant_index,
                     -1};
}

/// Make an OccPosition for an atom in a molecule
OccPosition OccSystem::make_atom_position(
    xtal::UnitCellCoord const &integral_site_coordinate,
    std::string orientation_name, Index atom_position_index) const {
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_position: Invalid "
        "integral_site_coordinate");
  }

  std::vector<xtal::Molecule> const &occupant_dof =
      this->prim->basis()[b].occupant_dof();
  Index occupant_index =
      find_index(this->prim->unique_names()[b], orientation_name);
  if (occupant_index < 0 || occupant_index >= occupant_dof.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_position: Invalid "
        "orientation_name");
  }

  xtal::Molecule const &mol = occupant_dof[occupant_index];
  if (atom_position_index < 0 || atom_position_index >= mol.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atomic_position: Invalid "
        "atom_position_index");
  }

  return OccPosition{false, true, integral_site_coordinate, occupant_index,
                     atom_position_index};
}

OccPosition OccSystem::make_molecule_in_reservoir_position(
    Index chemical_index) const {
  if (chemical_index < 0 || chemical_index >= this->chemical_name_list.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_in_reservoir_position: Invalid "
        "chemical_index");
  }
  return OccPosition{true, false, xtal::UnitCellCoord{0, 0, 0, 0},
                     chemical_index, -1};
}

/// Make an OccPosition that indicates occupant in the reservoir
OccPosition OccSystem::make_molecule_in_reservoir_position(
    std::string chemical_name) const {
  Index chemical_index = find_index(this->chemical_name_list, chemical_name);
  if (chemical_index < 0 || chemical_index >= this->chemical_name_list.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_molecule_in_reservoir_position: Invalid "
        "chemical_name");
  }
  return OccPosition{true, false, xtal::UnitCellCoord{0, 0, 0, 0},
                     chemical_index, -1};
}

/// Make an OccPosition that indicates atomic molecule in the reservoir
OccPosition OccSystem::make_atom_in_reservoir_position(
    std::string chemical_name) const {
  Index occupant_index = find_index(this->chemical_name_list, chemical_name);
  if (occupant_index < 0 || occupant_index >= this->chemical_name_list.size()) {
    throw std::runtime_error(
        "Error in OccSystem::make_atom_in_reservoir_position: Invalid "
        "OccPosition "
        "chemical_name");
  }
  return OccPosition{true, true, xtal::UnitCellCoord{0, 0, 0, 0},
                     occupant_index, 0};
}

/// \brief Populate OccPosition vector based on cluster and occupation
///
/// \param positions, OccPosition vector to be populated.
/// \param count, Eigen::VectorXi used for checking chemical counts.
/// \param cluster, Cluster where event is occuring
/// \param occ_init, Occupant indices before event occurs
/// \param occ_final, Occupant indices after event occurs
/// \param require_atom_conservation, If true, atoms on the cluster
///     before and after the event are conserved. Current implementation
///     only allows `require_atom_conservation==true`.
///
/// Notes:
/// - `positions` is populated with an OccPosition for each atom in each
///   molecule on each site in `occ_init`.
/// - To be implemented: include reservoir OccPosition determined from
///   difference between `occ_init` and `occ_final`.
void OccSystem::make_occ_positions(std::vector<OccPosition> &positions,
                                   Eigen::VectorXi &count,
                                   clust::IntegralCluster const &cluster,
                                   std::vector<int> const &occ_init,
                                   std::vector<int> const &occ_final,
                                   bool require_atom_conservation) const {
  // initial positions
  positions.clear();
  Index i = 0;
  for (auto const &integral_site_coordinate : cluster) {
    Index b = integral_site_coordinate.sublattice();
    xtal::Site site = prim->basis()[b];
    xtal::Molecule const &mol = site.occupant_dof()[occ_init[i]];
    for (Index j = 0; j < mol.size(); ++j) {
      positions.emplace_back(false, true, integral_site_coordinate, occ_init[i],
                             j);
    }
    ++i;
  }

  if (!require_atom_conservation) {
    throw std::runtime_error(
        "Error in OccSystem::make_occ_positions: non-conserved occ_position is "
        "not "
        "yet implemented.");
    // this will use parameter `count`
    // if non-conserved indivisible molecules, add mol-in-reservoir
    // add remaining as atom-in-reservoir
  }
}

Eigen::Vector3d OccSystem::get_cartesian_coordinate(
    OccPosition const &occ_position) const {
  auto const &pos = occ_position;
  if (pos.is_in_reservoir) {
    throw std::runtime_error(
        "Error in OccSystem::get_cartesian_coordinate: OccPosition is in "
        "reservoir");
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
  Eigen::Vector3d coord =
      pos.integral_site_coordinate.coordinate(*this->prim).cart();
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

/// \brief Return reference to molecule occupant
///
/// Note:
/// - For reservoir positions, this returns the first xtal::Molecule orientation
///   with matching chemical name
xtal::Molecule const &OccSystem::get_occupant(OccPosition const &pos) const {
  if (pos.is_in_reservoir) {
    std::string chemical_name = this->get_chemical_name(pos);
    for (Index b = 0; b < prim->basis().size(); ++b) {
      for (auto const &mol : prim->basis()[b].occupant_dof()) {
        if (mol.name() == chemical_name) {
          return mol;
        }
      }
    }
    throw std::runtime_error(
        "Error in OccSystem::get_occupant: Invalid OccPosition "
        "chemical_name_list");
  }

  Index b = pos.integral_site_coordinate.sublattice();
  if (b < 0 || b >= this->prim->basis().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_occupant: Invalid OccPosition "
        "sublattice");
  }
  xtal::Site const &site = this->prim->basis()[b];
  if (pos.occupant_index < 0 ||
      pos.occupant_index >= site.occupant_dof().size()) {
    throw std::runtime_error(
        "Error in OccSystem::get_occupant: Invalid OccPosition "
        "occupant_index");
  }
  return site.occupant_dof()[pos.occupant_index];
}

/// \brief Count molecules occupying cluster sites
///
/// \param count, Eigen::VectorXi used for checking atom count.
///     Final value will be store here, with entries in order of
///     this->chemical_name_list. Will be resized if necessary.
/// \param cluster, Cluster where event is occuring
/// \param occ, Occupant indices on each cluster site
///
/// Note:
/// - This will segfault if parameters are not sized correctly
///   or give out of range values correctly
void OccSystem::molecule_count(Eigen::VectorXi &count,
                               clust::IntegralCluster const &cluster,
                               std::vector<int> const &occ) const {
  count.resize(chemical_name_list.size());
  count.setZero();
  Index b;
  Index i = 0;
  for (auto const &site : cluster) {
    b = site.sublattice();
    count(occupant_to_chemical_index[b][occ[i]]) += 1;
    ++i;
  }
  return;
}

/// \brief Count molecules (by orientation) occupying cluster sites
///
/// \param count, Eigen::VectorXi used for checking atom count.
///     Final value will be store here, with entries in order of
///     this->chemical_name_list. Will be resized if necessary.
/// \param cluster, Cluster where event is occuring
/// \param occ, Occupant indices on each cluster site
///
/// Note:
/// - This will segfault if parameters are not sized correctly
///   or give out of range values correctly
void OccSystem::orientation_count(Eigen::VectorXi &count,
                                  clust::IntegralCluster const &cluster,
                                  std::vector<int> const &occ) const {
  count.resize(orientation_name_list.size());
  count.setZero();
  Index b;
  Index i = 0;
  for (auto const &site : cluster) {
    b = site.sublattice();
    count(occupant_to_orientation_index[b][occ[i]]) += 1;
    ++i;
  }
  return;
}

/// \brief Count atoms occupying cluster sites
///
/// \param count, Eigen::VectorXi used for checking atom count.
///     Final value will be stored here, with entries in order
///     of this->atom_name_list. Will be resized if necessary.
/// \param cluster, Cluster where event is occuring
/// \param occ, Occupant indices on each cluster site
///
/// Note:
/// - This will segfault if parameters are not sized correctly
///   or give out of range values correctly
void OccSystem::atom_count(Eigen::VectorXi &count,
                           clust::IntegralCluster const &cluster,
                           std::vector<int> const &occ) const {
  count.resize(atom_name_list.size());
  count.setZero();
  Index b;
  Index i = 0;
  for (auto const &site : cluster) {
    b = site.sublattice();
    for (auto const &atom_name_index : atom_position_to_name_index[b][occ[i]]) {
      count(atom_name_index) += 1;
    }
    ++i;
  }
  return;
}

/// \brief Check if potential event is conserves molecule chemical type
///
/// \param count, Eigen::VectorXi used for checking chemical count.
///     Final value will be set to count_final - count_init, with
///     entries in order of this->chemical_name_list.
/// \param cluster, Cluster where event is occuring
/// \param occ_init, Occupant indices before event occurs
/// \param occ_final, Occupant indices after event occurs
///
/// Note:
/// - This will segfault if parameters are not sized correctly
///   or give out of range values correctly
bool OccSystem::is_molecule_conserving(
    Eigen::VectorXi &count, clust::IntegralCluster const &cluster,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final) const {
  count.resize(chemical_name_list.size());
  count.setZero();
  Index b;
  Index i = 0;
  for (auto const &site : cluster) {
    b = site.sublattice();
    count(occupant_to_chemical_index[b][occ_final[i]]) += 1;
    count(occupant_to_chemical_index[b][occ_init[i]]) -= 1;
    ++i;
  }
  return !count.any();
}

/// \brief Check if potential event is conserves atom chemical type
///
/// \param count, Eigen::VectorXi used for checking atom count.
///     Final value will be set to count_final - count_init, with
///     entries in order of this->atom_name_list.
/// \param cluster, Cluster where event is occuring
/// \param occ_init, Occupant indices before event occurs
/// \param occ_final, Occupant indices after event occurs
///
/// Note:
/// - This will segfault if parameters are not sized correctly
///   or give out of range values correctly
bool OccSystem::is_atom_conserving(Eigen::VectorXi &count,
                                   clust::IntegralCluster const &cluster,
                                   std::vector<int> const &occ_init,
                                   std::vector<int> const &occ_final) const {
  count.resize(atom_name_list.size());
  count.setZero();
  Index b;
  Index i = 0;
  for (auto const &site : cluster) {
    b = site.sublattice();
    for (auto const &atom_name_index :
         atom_position_to_name_index[b][occ_final[i]]) {
      count(atom_name_index) += 1;
    }
    for (auto const &atom_name_index :
         atom_position_to_name_index[b][occ_init[i]]) {
      count(atom_name_index) -= 1;
    }
    ++i;
  }
  return !count.any();
}

/// \brief Compares atom_name (for atoms) and/or chemical name (for molecule)
bool OccSystem::is_same_chemical_type(OccPosition const &pos1,
                                      OccPosition const &pos2) const {
  if (pos1.is_atom) {
    if (pos2.is_atom) {
      return get_atom_name_index(pos1) == get_atom_name_index(pos2);
    }
    return get_atom_name(pos1) == get_chemical_name(pos2);
  } else {
    if (pos2.is_atom) {
      return get_chemical_name(pos1) == get_atom_name(pos2);
    }
    return get_chemical_index(pos1) == get_chemical_index(pos2);
  }
}

/// \brief Check if each before and after position is the same chemical type
bool OccSystem::is_chemical_type_conserving(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  for (int i = 0; i < position_before.size(); ++i) {
    if (!is_same_chemical_type(position_before[i], position_after[i])) {
      return false;
    }
  }
  return true;
}

/// \brief Check for trajectories in which two atoms
///     directly exchange sites. Does not include atom-vacancy exchange.
bool OccSystem::is_direct_exchange(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  auto const &before = position_before;
  auto const &after = position_after;

  for (auto const &pos : before) {
    if (pos.is_in_reservoir || this->is_vacancy(pos)) {
      return false;
    }
  }
  for (auto const &pos : after) {
    if (pos.is_in_reservoir || this->is_vacancy(pos)) {
      return false;
    }
  }

  for (Index i = 0; i < before.size() - 1; ++i) {
    for (Index j = i + 1; j < after.size(); ++j) {
      if (before[i].integral_site_coordinate ==
              after[j].integral_site_coordinate &&
          before[i].atom_position_index == after[j].atom_position_index &&
          before[j].integral_site_coordinate ==
              after[i].integral_site_coordinate &&
          before[j].atom_position_index == after[i].atom_position_index) {
        return true;
      }
    }
  }
  return false;
}

/// \brief Return true if any site is a vacancy before and after the event
bool OccSystem::is_any_unchanging_vacant_site(
    clust::IntegralCluster const &cluster, std::vector<int> const &occ_init,
    std::vector<int> const &occ_final) const {
  Index i = 0;
  for (auto const &site : cluster) {
    Index b = site.sublattice();
    if (is_vacancy_list[occupant_to_chemical_index[b][occ_init[i]]] &&
        is_vacancy_list[occupant_to_chemical_index[b][occ_final[i]]]) {
      return true;
    }
    ++i;
  }
  return false;
}

/// \brief Return true if before and after position are both in reservoir
///
bool OccSystem::is_any_unchanging_reservoir_type(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  // check if before and after are in reservoir or before and after are same
  // site
  for (Index i = 0; i < position_before.size(); ++i) {
    if (position_before[i].is_in_reservoir &&
        position_after[i].is_in_reservoir) {
      return true;
    }
  }
  return false;
}

/// \brief Return true if, for any molecule, all atoms remain fixed
bool OccSystem::is_any_unchanging_molecule(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  // initialize site:<has_a_change> to false for all sites
  std::map<xtal::UnitCellCoord, bool> has_a_change;
  for (auto const &p1 : position_before) {
    if (p1.is_in_reservoir == false) {
      has_a_change[p1.integral_site_coordinate] = false;
    }
  }

  // iterate over atom trajectories and check if site changes
  // or occupant changes or position changes
  auto p2_it = position_after.begin();
  for (auto const &p1 : position_before) {
    auto const &p2 = *p2_it;
    if (p1.is_in_reservoir && p2.is_in_reservoir) {
      // pass
    } else if (p1.is_in_reservoir && !p2.is_in_reservoir) {
      has_a_change[p2.integral_site_coordinate] = true;
    } else if (!p1.is_in_reservoir && p2.is_in_reservoir) {
      has_a_change[p1.integral_site_coordinate] = true;
    } else if (p1.integral_site_coordinate != p2.integral_site_coordinate) {
      has_a_change[p1.integral_site_coordinate] = true;
      has_a_change[p2.integral_site_coordinate] = true;
    } else if (p1.occupant_index != p2.occupant_index) {
      has_a_change[p1.integral_site_coordinate] = true;
    } else if (p1.atom_position_index != p2.atom_position_index) {
      has_a_change[p1.integral_site_coordinate] = true;
    }

    ++p2_it;
  }

  for (auto const &pair : has_a_change) {
    if (pair.second == false) {
      return true;
    }
  }
  return false;
}

/// \brief Return true if, for any molecule, component atoms
///     begin/end in different molecules
bool OccSystem::is_molecule_breakup(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  std::map<xtal::UnitCellCoord, std::set<xtal::UnitCellCoord>> init_to_final;
  for (auto const &p1 : position_before) {
    if (p1.is_in_reservoir) {
      continue;
    }
    init_to_final[p1.integral_site_coordinate].clear();
  }

  for (auto const &p2 : position_after) {
    if (p2.is_in_reservoir) {
      continue;
    }
    init_to_final[p2.integral_site_coordinate].clear();
  }

  // iterate over trajectories and insert init site -> final site
  auto p2_it = position_after.begin();
  for (auto const &p1 : position_before) {
    auto const &p2 = *p2_it;
    if (p1.is_in_reservoir) {
      continue;
    }

    init_to_final[p1.integral_site_coordinate].insert(
        p2.integral_site_coordinate);
    ++p2_it;
  }

  // iterate over trajectories and insert init site <- final site
  auto p1_it = position_before.begin();
  for (auto const &p2 : position_after) {
    auto const &p1 = *p1_it;
    if (p2.is_in_reservoir) {
      continue;
    }

    init_to_final[p2.integral_site_coordinate].insert(
        p1.integral_site_coordinate);
    ++p1_it;
  }

  for (auto const &pair : init_to_final) {
    if (pair.second.size() > 1) {
      return true;
    }
  }
  return false;
}

/// \brief Return true if, for any indivisible molecule, component atoms
///     begin/end in different molecules
bool OccSystem::is_indivisible_molecule_breakup(
    std::vector<OccPosition> const &position_before,
    std::vector<OccPosition> const &position_after) const {
  std::map<xtal::UnitCellCoord, std::set<xtal::UnitCellCoord>> init_to_final;

  auto _reset = [&](std::vector<OccPosition> const &position_before,
                    std::vector<OccPosition> const &position_after) {
    init_to_final.clear();
    for (auto const &p1 : position_before) {
      if (p1.is_in_reservoir || !is_indivisible(p1)) {
        continue;
      }
      init_to_final[p1.integral_site_coordinate].clear();
    }
  };

  auto _count = [&](std::vector<OccPosition> const &position_before,
                    std::vector<OccPosition> const &position_after) {
    auto p2_it = position_after.begin();
    for (auto const &p1 : position_before) {
      auto const &p2 = *p2_it;
      if (p1.is_in_reservoir || !is_indivisible(p1)) {
        continue;
      }

      init_to_final[p1.integral_site_coordinate].insert(
          p2.integral_site_coordinate);
      ++p2_it;
    }
  };

  auto _check = [&]() {
    for (auto const &pair : init_to_final) {
      if (pair.second.size() > 1) {
        return true;
      }
    }
    return false;
  };

  // -- check forward --
  _reset(position_before, position_after);
  _count(position_before, position_after);
  if (_check()) {
    return true;
  }

  // -- check reverse --
  _reset(position_after, position_before);
  _count(position_after, position_before);
  return _check();
}

/// \brief Check if a molecule is contained in a list, in given orientation
bool is_contained_in_this_orientation(
    std::vector<xtal::Molecule> const &molecule_list,
    xtal::Molecule const &molecule, double tol) {
  if (molecule.name().empty()) {
    throw std::runtime_error("Error: molecule has empty name");
  }
  auto it = std::find_if(
      molecule_list.begin(), molecule_list.end(),
      [&](xtal::Molecule const &x) { return x.identical(molecule, tol); });
  if (it != molecule_list.end()) {
    if (it->name() != molecule.name()) {
      throw std::runtime_error(
          "Error: equivalent molecules have different names");
    }
    return true;
  }
  return false;
}

/// \brief Check if a molecule is contained in a list, in any orientation
bool is_contained_in_any_orientation(
    std::vector<xtal::Molecule> const &molecule_list,
    xtal::Molecule const &molecule,
    std::vector<xtal::SymOp> const &factor_group, double tol) {
  for (auto const &op : factor_group) {
    xtal::Molecule tmol = sym::copy_apply(op, molecule);
    if (is_contained_in_this_orientation(molecule_list, tmol, tol)) {
      return true;
    }
  }
  return false;
}

/// \brief Check that all molecules have non-empty names and that
///     symmetrically equivalent molecules have the same name
bool is_valid_molecule_naming(xtal::BasicStructure const &prim,
                              std::vector<xtal::SymOp> const &factor_group) {
  try {
    auto molecule_list = molecule_list_single_orientation(prim, factor_group);
  } catch (std::exception &e) {
    return false;
  }
  return true;
}

/// \brief Generate a list of all molecule orientations in a prim
std::vector<xtal::Molecule> molecule_list_all_orientations(
    xtal::BasicStructure const &prim) {
  double tol = prim.lattice().tol();
  std::vector<xtal::Molecule> molecule_list;
  for (auto const &site : prim.basis()) {
    for (auto const &mol : site.occupant_dof()) {
      if (!is_contained_in_this_orientation(molecule_list, mol, tol)) {
        molecule_list.emplace_back(mol);
      }
    }
  }
  return molecule_list;
}

/// \brief Generate a list of symmetrically unique molecules in a prim
std::vector<xtal::Molecule> molecule_list_single_orientation(
    xtal::BasicStructure const &prim,
    std::vector<xtal::SymOp> const &factor_group) {
  double tol = prim.lattice().tol();
  std::vector<xtal::Molecule> molecule_list;
  for (auto const &site : prim.basis()) {
    for (auto const &mol : site.occupant_dof()) {
      if (!is_contained_in_any_orientation(molecule_list, mol, factor_group,
                                           tol)) {
        molecule_list.emplace_back(mol);
      }
    }
  }
  return molecule_list;
}

/// \brief Generate a list of names of molecules in a prim, using
///     the `prim->unique_names()`
std::vector<std::string> make_orientation_name_list(
    xtal::BasicStructure const &prim) {
  std::vector<xtal::Molecule> molecule_list =
      molecule_list_all_orientations(prim);
  std::vector<std::string> orientation_name_list;
  for (auto const &mol : molecule_list) {
    orientation_name_list.push_back(orientation_name(mol, prim));
  }
  return orientation_name_list;
}

/// \brief Generate a list of names of symmetrically unique molecules
///     from a list of all molecule orientations
std::vector<std::string> make_chemical_name_list(
    xtal::BasicStructure const &prim,
    std::vector<xtal::SymOp> const &factor_group) {
  std::vector<xtal::Molecule> molecule_list =
      molecule_list_single_orientation(prim, factor_group);
  std::vector<std::string> chemical_name_list;
  for (auto const &mol : molecule_list) {
    chemical_name_list.push_back(mol.name());
  }
  return chemical_name_list;
}

/// \brief Generate a list of unique atom names (sorted)
std::vector<std::string> make_atom_name_list(xtal::BasicStructure const &prim) {
  std::vector<std::string> atom_name_list;

  Index b = 0;
  for (auto const &site : prim.basis()) {
    Index occupant_index = 0;
    for (auto const &mol : site.occupant_dof()) {
      for (auto const &atom : mol.atoms()) {
        Index atom_name_index = find_index(atom_name_list, atom.name());
        if (atom_name_index == atom_name_list.size()) {
          atom_name_index = atom_name_list.size();
          atom_name_list.push_back(atom.name());
        }
      }
      ++occupant_index;
    }

    ++b;
  }
  std::sort(atom_name_list.begin(), atom_name_list.end());
  return atom_name_list;
}

/// \brief Generate a list of molecules in a prim, exclude symmetrically
///     equivalent orientations
std::vector<xtal::Molecule> molecule_list_single_orientation(
    std::vector<xtal::Molecule> const &molecule_list_all_orientations,
    std::vector<xtal::SymOp> const &factor_group, double tol) {
  std::vector<xtal::Molecule> molecule_list;
  for (auto const &mol : molecule_list_all_orientations) {
    if (!is_contained_in_any_orientation(molecule_list, mol, factor_group,
                                         tol)) {
      molecule_list.emplace_back(mol);
    }
  }
  return molecule_list;
}

/// \brief Get the unique name of a specific molecule orientation
std::string orientation_name(xtal::Molecule const &molecule,
                             xtal::BasicStructure const &prim) {
  double tol = prim.lattice().tol();
  if (prim.unique_names().size() != prim.basis().size()) {
    throw std::runtime_error("Error in orientation_name: basis size mismatch");
  }

  auto basis_unique_names_it = prim.unique_names().begin();
  for (auto const &site : prim.basis()) {
    if (basis_unique_names_it->size() != site.occupant_dof().size()) {
      throw std::runtime_error(
          "Error in orientation_name: occupant size mismatch");
    }
    auto occupant_dof_unique_names_it = basis_unique_names_it->begin();
    for (auto const &mol : site.occupant_dof()) {
      if (molecule.identical(mol, tol)) {
        return *occupant_dof_unique_names_it;
      }
      ++occupant_dof_unique_names_it;
    }
    ++basis_unique_names_it;
  }
  throw std::runtime_error(
      "Error in orientation_name: molecule not found in prim");
}

}  // namespace occ_events
}  // namespace CASM
