#ifndef CASM_occ_events_OccSystem
#define CASM_occ_events_OccSystem

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class Molecule;
struct SymOp;
class UnitCellCoord;
}  // namespace xtal

namespace clust {
class IntegralCluster;
}  // namespace clust

namespace occ_events {

struct OccPosition;

/// \brief Defines the system for OccPosition / OccTrajectory / OccEvent
///     and provides lookup tables and methods for determining properties
///
/// Notes:
/// - Symmetry is only necessary to determine `chemical_name_list`, which
///   can be done using the `make_chemical_name_list` function. From there,
///   it is presumed that `xtal::Molecule::name()` is the same for all
///   `xtal::Molecule` that have the same chemical type.
/// - The `orientation_name_list` names all non-identical xtal::Molecule.
///   These names match `xtal::BasicStructure::unique_names()`.
/// - The `atom_name_list` names all atom components of xtal::Molecule.
///   These names match `xtal::AtomPosition::name()`.
struct OccSystem {
  OccSystem(std::shared_ptr<xtal::BasicStructure const> const &_prim,
            std::vector<std::string> const &_chemical_name_list,
            std::set<std::string> const &_vacancy_name_list =
                std::set<std::string>({"Va", "VA", "va"}));

  /// \brief The prim
  std::shared_ptr<xtal::BasicStructure const> prim;

  /// \brief Names of the unique chemical components
  ///
  /// Notes:
  /// - When `is_in_reservoir==true` OccPosition::occupant_index
  ///   is an index into this / chemical_list
  /// - These come from Molecule::name
  /// - `chemical_index` indicates an index into this list
  std::vector<std::string> chemical_name_list;

  /// \brief True if Molecule is indivisible
  ///
  /// Notes:
  /// - Use `chemical_index` to index into this list
  std::vector<bool> is_indivisible_chemical_list;

  /// \brief Chemical names that should be counted as vacancies
  std::set<std::string> vacancy_name_list;

  /// \brief Check if molecule is vacancy (via chemical_index)
  std::vector<bool> is_vacancy_list;

  /// \brief Names of the unique atomic components (by name)
  ///
  /// Notes:
  /// - `atom_name_index` indicates an index into this list
  /// - This only distinguishes by atom `name`, not atom
  ///   properties.
  /// - TODO: `atom_type_list` to distinguish by fixed properties
  std::vector<std::string> atom_name_list;

  /// \brief Names of the unique molecular orientations
  ///
  /// Notes:
  /// - These come from `prim->unique_names()`
  /// - `orientation_index` indicates an index into this list
  std::vector<std::string> orientation_name_list;

  /// \brief Conversion table of sublattice & occupant index to chemical index
  ///
  /// Usage:
  /// - atom_position_to_name_index[b][occupant_index][atom_position_index] ->
  /// atom_name_index
  std::vector<std::vector<std::vector<int>>> atom_position_to_name_index;

  /// \brief Conversion table of sublattice & occupant index to chemical index
  ///
  /// Usage:
  /// - occupant_to_chemical_index[b][occupant_index] -> chemical_index
  std::vector<std::vector<int>> occupant_to_chemical_index;

  /// \brief Conversion table of sublattice & occupant index to orientation
  /// index
  ///
  /// Usage:
  /// - occupant_to_orientation_index[b][occupant_index] -> orientation_index
  /// - -1 if invalid
  std::vector<std::vector<int>> occupant_to_orientation_index;

  /// \brief Conversion table of sublattice & orientation index to occupant
  /// index
  ///
  /// Usage:
  /// - occupant_to_orientation_index[b][orientation_index] -> occupant_index
  /// - -1 if invalid
  std::vector<std::vector<int>> orientation_to_occupant_index;

  // --- OccPosition factory functions ---

  OccPosition make_molecule_position(
      xtal::UnitCellCoord const &integral_site_coordinate,
      Index occupant_index) const;

  OccPosition make_atom_position(
      xtal::UnitCellCoord const &integral_site_coordinate, Index occupant_index,
      Index atom_position_index = 0) const;

  OccPosition make_molecule_position(
      xtal::UnitCellCoord const &integral_site_coordinate,
      std::string orientation_name) const;

  OccPosition make_atom_position(
      xtal::UnitCellCoord const &integral_site_coordinate,
      std::string orientation_name, Index atom_position_index = 0) const;

  OccPosition make_molecule_in_reservoir_position(Index chemical_index) const;

  OccPosition make_molecule_in_reservoir_position(
      std::string chemical_name) const;

  OccPosition make_atom_in_reservoir_position(std::string chemical_name) const;

  void make_occ_positions(std::vector<OccPosition> &positions,
                          Eigen::VectorXi &count,
                          clust::IntegralCluster const &cluster,
                          std::vector<int> const &occ_init,
                          std::vector<int> const &occ_final,
                          bool require_atom_conservation) const;

  // --- Get OccPosition info ---

  std::string get_chemical_name(
      xtal::UnitCellCoord const &integral_site_coordinate,
      Index occupant_index) const {
    return chemical_name_list
        [occupant_to_chemical_index[integral_site_coordinate.sublattice()]
                                   [occupant_index]];
  }

  std::string get_orientation_name(
      xtal::UnitCellCoord const &integral_site_coordinate,
      Index occupant_index) const {
    return orientation_name_list
        [occupant_to_orientation_index[integral_site_coordinate.sublattice()]
                                      [occupant_index]];
  }

  std::string get_atom_name(xtal::UnitCellCoord const &integral_site_coordinate,
                            Index occupant_index,
                            Index atom_position_index) const {
    return atom_name_list
        [atom_position_to_name_index[integral_site_coordinate.sublattice()]
                                    [occupant_index][atom_position_index]];
  }

  std::vector<std::string> get_atom_names(
      xtal::UnitCellCoord const &integral_site_coordinate,
      Index occupant_index) const {
    std::vector<std::string> atom_names;
    auto const &indices =
        atom_position_to_name_index[integral_site_coordinate.sublattice()]
                                   [occupant_index];
    for (auto const &atom_name_index : indices) {
      atom_names.push_back(atom_name_list[atom_name_index]);
    }
    return atom_names;
  }

  /// Valid if p is valid
  std::string get_chemical_name(OccPosition const &p) const {
    return chemical_name_list[get_chemical_index(p)];
  }

  /// Valid if p.is_in_reservoir==false
  std::string get_orientation_name(OccPosition const &p) const {
    return orientation_name_list[get_orientation_index(p)];
  }

  /// Valid if p.is_atom==true
  std::string get_atom_name(OccPosition const &p) const {
    return atom_name_list[get_atom_name_index(p)];
  }

  /// Valid if p.is_in_reservoir==false
  Index get_sublattice_index(OccPosition const &p) const {
    return p.integral_site_coordinate.sublattice();
  }

  /// Valid if p is valid
  Index get_chemical_index(OccPosition const &p) const {
    if (p.is_in_reservoir) {
      return p.occupant_index;
    }
    return occupant_to_chemical_index[get_sublattice_index(p)]
                                     [p.occupant_index];
  }

  /// Valid if p.is_in_reservoir==false
  Index get_orientation_index(OccPosition const &p) const {
    return occupant_to_orientation_index[get_sublattice_index(p)]
                                        [p.occupant_index];
  }

  /// Valid if p.is_atom==true
  Index get_atom_name_index(OccPosition const &p) const {
    return atom_position_to_name_index[get_sublattice_index(p)]
                                      [p.occupant_index][p.atom_position_index];
  }

  /// Valid if p.is_in_reservoir==false
  Eigen::Vector3d get_cartesian_coordinate(OccPosition const &p) const;

  /// Return true if molecule is indivisible
  bool is_indivisible(OccPosition const &p) const {
    return is_indivisible_chemical_list[get_chemical_index(p)];
  }

  bool is_vacancy(OccPosition const &p) const {
    return is_vacancy_list[get_chemical_index(p)];
  }

  /// \brief Return reference to molecule occupant
  xtal::Molecule const &get_occupant(OccPosition const &p) const;

  // --- Occupation checks ---

  /// \brief Count molecules occupying cluster sites
  void molecule_count(Eigen::VectorXi &count,
                      clust::IntegralCluster const &cluster,
                      std::vector<int> const &occ) const;

  /// \brief Count molecules (by orientation) occupying cluster sites
  void orientation_count(Eigen::VectorXi &count,
                         clust::IntegralCluster const &cluster,
                         std::vector<int> const &occ) const;

  /// \brief Count atoms occupying cluster sites
  void atom_count(Eigen::VectorXi &count, clust::IntegralCluster const &cluster,
                  std::vector<int> const &occ) const;

  /// \brief Check if potential event is conserves molecule chemical type
  bool is_molecule_conserving(Eigen::VectorXi &count,
                              clust::IntegralCluster const &cluster,
                              std::vector<int> const &occ_init,
                              std::vector<int> const &occ_final) const;

  /// \brief Check if potential event is conserves atom chemical type
  bool is_atom_conserving(Eigen::VectorXi &count,
                          clust::IntegralCluster const &cluster,
                          std::vector<int> const &occ_init,
                          std::vector<int> const &occ_final) const;

  /// \brief Return true if any site is a vacancy before and after the event
  bool is_any_unchanging_vacant_site(clust::IntegralCluster const &cluster,
                                     std::vector<int> const &occ_init,
                                     std::vector<int> const &occ_final) const;

  // --- Trajectory checks ---

  /// \brief Compares atom_name (for atoms) and/or chemical name (for molecule)
  bool is_same_chemical_type(OccPosition const &pos1,
                             OccPosition const &pos2) const;

  /// \brief Check if each before and after position is the same chemical type
  bool is_chemical_type_conserving(
      std::vector<OccPosition> const &position_before,
      std::vector<OccPosition> const &position_after) const;

  /// \brief Check for trajectories in which two atoms
  ///     directly exchange sites. Does not include atom-vacancy exchange.
  bool is_direct_exchange(std::vector<OccPosition> const &position_before,
                          std::vector<OccPosition> const &position_after) const;

  /// \brief Return true if before and after position are both in reservoir
  bool is_any_unchanging_reservoir_type(
      std::vector<OccPosition> const &position_before,
      std::vector<OccPosition> const &position_after) const;

  /// \brief Return true if, for any molecule, all atoms remain fixed
  bool is_any_unchanging_molecule(
      std::vector<OccPosition> const &position_before,
      std::vector<OccPosition> const &position_after) const;

  /// \brief Return true if, for any molecule, component atoms
  ///     begin/end in different molecules
  bool is_molecule_breakup(
      std::vector<OccPosition> const &position_before,
      std::vector<OccPosition> const &position_after) const;

  /// \brief Return true if, for any indivisible molecule, component atoms
  ///     begin/end in different molecules
  bool is_indivisible_molecule_breakup(
      std::vector<OccPosition> const &position_before,
      std::vector<OccPosition> const &position_after) const;
};

/// \brief Check if a molecule is contained in a list, in given orientation
bool is_contained_in_this_orientation(
    std::vector<xtal::Molecule> const &molecule_list,
    xtal::Molecule const &molecule, double tol);

/// \brief Check if a molecule is contained in a list, in any orientation
bool is_contained_in_any_orientation(
    std::vector<xtal::Molecule> const &molecule_list,
    xtal::Molecule const &molecule,
    std::vector<xtal::SymOp> const &factor_group, double tol);

/// \brief Check that all molecules have non-empty names and that
///     symmetrically equivalent molecules have the same name
bool is_valid_molecule_naming(xtal::BasicStructure const &prim,
                              std::vector<xtal::SymOp> const &factor_group);

/// \brief Generate a list of all molecule orientations in a prim
std::vector<xtal::Molecule> molecule_list_all_orientations(
    xtal::BasicStructure const &prim);

/// \brief Generate a list of molecules in a prim, exclude symmetrically
///     equivalent orientations
std::vector<xtal::Molecule> molecule_list_single_orientation(
    xtal::BasicStructure const &prim,
    std::vector<xtal::SymOp> const &factor_group);

/// \brief Generate a list of symmetrically unique molecules from a
///     list of all molecule orientations
std::vector<xtal::Molecule> molecule_list_single_orientation(
    std::vector<xtal::Molecule> const &molecule_list_all_orientations,
    std::vector<xtal::SymOp> const &factor_group, double tol);

/// \brief Generate a list of names of molecules in a prim, using
///     the `prim->unique_names()`
std::vector<std::string> make_orientation_name_list(
    xtal::BasicStructure const &prim);

/// \brief Generate a list of names of symmetrically unique molecules
///     from a list of all molecule orientations
std::vector<std::string> make_chemical_name_list(
    xtal::BasicStructure const &prim,
    std::vector<xtal::SymOp> const &factor_group);

/// \brief Generate a list of unique atom names
std::vector<std::string> make_atom_name_list(xtal::BasicStructure const &prim);

/// \brief Get the unique name of a specific molecule orientation
std::string orientation_name(xtal::Molecule const &molecule,
                             xtal::BasicStructure const &prim);

}  // namespace occ_events
}  // namespace CASM

#endif
