#ifndef CASM_occ_events_OccPosition
#define CASM_occ_events_OccPosition

#include "casm/configuration/occ_events/definitions.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace occ_events {

/// An atom or molecule position
///
/// - Allows specifying a single-atom or multi-atom molecule
/// - Allows specifying a single atom inside a multi-atom molecule
/// - Allows specifying atoms/molecules "in the resevoir" for
///   describing (semi) grand canonical Monte Carlo events
struct OccPosition
    : public Comparisons<xtal::Translatable<CRTPBase<OccPosition>>> {
  OccPosition(bool _is_in_resevoir, bool _is_atom,
              xtal::UnitCellCoord const &_integral_site_coordinate,
              Index _occupant_index, Index _atom_position_index);

  static OccPosition molecule_in_resevoir(Index _occupant_index);

  /// \brief Only used for single-atom molecules
  static OccPosition atom_in_resevoir(Index _occupant_index);

  static OccPosition molecule(
      xtal::UnitCellCoord const &_integral_site_coordinate,
      Index _occupant_index);

  static OccPosition atom(xtal::UnitCellCoord const &_integral_site_coordinate,
                          Index _occupant_index, Index _atom_position_index);

  /// If true, indicates molecule/atom in resevoir. If false,
  /// indicates a molecule/atom on integral_site_coordinate
  bool is_in_resevoir;

  /// If true, indicates this tracks an atom position. If false, then
  /// this tracks a molecule position.
  bool is_atom;

  /// - If is_in_resevoir==false: Integral coordinates of site containing
  /// occupant
  /// - Otherwise: Ignored
  xtal::UnitCellCoord integral_site_coordinate;

  /// - If is_in_resevoir==false: Index of Molecule in Site::occupant_dof()
  ///   for sublattice specified by `integral_site_coordinate`
  /// - If is_in_resevoir==true: Index of molecule name in
  ///   `OccSystem::resevoir_components`.
  Index occupant_index;

  /// - If is_atom==true && is_in_resevoir==false: Index of atom position
  ///   in Molecule::atoms()
  /// - If is_atom==true && is_in_resevoir==true: When indicating a Molecule
  ///   used to represent a single atom, set to zero 0.
  Index atom_position_index;

  /// Translate theOccPosition by a UnitCell translation
  OccPosition &operator+=(xtal::UnitCell trans);

  /// Compare (integral_site_coordinate, occupant_index, atom_position_index)
  bool operator<(OccPosition const &B) const;
};

/// \brief Apply SymOp to OccPosition
OccPosition &apply(OccEventRep const &rep, OccPosition &occ_position);

/// \brief Apply SymOp to OccPosition
OccPosition copy_apply(OccEventRep const &rep, OccPosition occ_position);

}  // namespace occ_events
}  // namespace CASM

#endif
