#include "casm/configuration/UnitCellCoordRep.hh"

#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace config {

// --- Inline definitions ---

/// \brief Apply symmetry to UnitCellCoord
///
/// \param rep UnitCellCoordRep representation of the symmetry operation
/// \param integral_site_coordinate Coordinate being transformed
///
UnitCellCoord &apply(UnitCellCoordRep const &rep,
                     UnitCellCoord &integral_site_coordinate) {
  integral_site_coordinate = copy_apply(rep, integral_site_coordinate);
  return integral_site_coordinate;
}

/// \brief Apply symmetry to UnitCellCoord
///
/// \param rep UnitCellCoordRep representation of the symmetry operation
/// \param integral_site_coordinate Coordinate being transformed
///
UnitCellCoord copy_apply(UnitCellCoordRep const &rep,
                         UnitCellCoord integral_site_coordinate) {
  UnitCell unitcell_indices =
      rep.point_matrix * integral_site_coordinate.unitcell() +
      rep.unitcell_indices[integral_site_coordinate.sublattice()];
  Index sublattice_index =
      rep.sublattice_index[integral_site_coordinate.sublattice()];
  return UnitCellCoord(sublattice_index, unitcell_indices);
}

}  // namespace config
}  // namespace CASM
