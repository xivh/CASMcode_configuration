#ifndef CASM_config_Supercell
#define CASM_config_Supercell

#include "casm/configuration/Prim.hh"
#include "casm/configuration/SupercellSymInfo.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace config {

/// \brief Specifies all the structural and symmetry information common for all
/// configurations with the same supercell. All members are const.
struct Supercell : public Comparisons<CRTPBase<Supercell>> {
  Supercell(std::shared_ptr<Prim const> const &_prim,
            Lattice const &_superlattice,
            Index max_n_translation_permutations = 100);
  Supercell(std::shared_ptr<Prim const> const &_prim,
            Superlattice const &_superlattice,
            Index max_n_translation_permutations = 100);
  Supercell(std::shared_ptr<Prim const> const &_prim,
            Eigen::Matrix3l const &_superlattice_matrix,
            Index max_n_translation_permutations = 100);

  /// \brief Species the primitive crystal structure (lattice and basis) and
  /// allowed degrees of freedom (DoF), and also symmetry representations
  /// used for all configurations with the same prim
  std::shared_ptr<Prim const> const prim;

  /// Couples the primitive lattice to the supercell lattice, and knows the
  /// transformation matrix
  Superlattice const superlattice;

  /// \brief Converts between ijk (UnitCell) values and their corresponding
  /// index in an unrolled vector
  xtal::UnitCellIndexConverter const unitcell_index_converter;

  /// \brief Converts between bijk (UnitCellCoord) values and their
  /// corresponding linear index
  xtal::UnitCellCoordIndexConverter const unitcellcoord_index_converter;

  /// \brief Holds symmetry representations used for all configurations with
  /// the same supercell
  SupercellSymInfo const sym_info;

  /// \brief Less than comparison of Supercell
  bool operator<(Supercell const &B) const;

 private:
  friend struct Comparisons<CRTPBase<Supercell>>;

  /// \brief Equality comparison of Supercell
  bool eq_impl(Supercell const &rhs) const;
};

struct CompareSharedSupercell {
  bool operator()(std::shared_ptr<Supercell const> const &lhs,
                  std::shared_ptr<Supercell const> const &rhs) const {
    return *lhs < *rhs;
  }
};

}  // namespace config
}  // namespace CASM

#endif
