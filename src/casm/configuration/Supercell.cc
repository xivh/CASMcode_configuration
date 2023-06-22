#include "casm/configuration/Supercell.hh"

#include "casm/configuration/SupercellSymInfo.hh"

namespace CASM {
namespace config {

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Lattice const &_superlattice,
                     Index max_n_translation_permutations)
    : Supercell(_prim,
                Superlattice(_prim->basicstructure->lattice(), _superlattice),
                max_n_translation_permutations) {}

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Superlattice const &_superlattice,
                     Index max_n_translation_permutations)
    : prim(_prim),
      superlattice(_superlattice),
      unitcell_index_converter(superlattice.transformation_matrix_to_super()),
      unitcellcoord_index_converter(
          superlattice.transformation_matrix_to_super(),
          prim->basicstructure->basis().size()),
      sym_info(prim, superlattice, unitcell_index_converter,
               unitcellcoord_index_converter, max_n_translation_permutations) {}

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Eigen::Matrix3l const &_superlattice_matrix,
                     Index max_n_translation_permutations)
    : Supercell(
          _prim,
          Superlattice(_prim->basicstructure->lattice(), _superlattice_matrix),
          max_n_translation_permutations) {}

/// \brief Less than comparison of Supercell
bool Supercell::operator<(Supercell const &B) const {
  if (prim != B.prim) {
    throw std::runtime_error(
        "Error using Supercell::operator<(Supercell const &B): "
        "Only Supercell with shared prim may be compared this way.");
  }
  if (superlattice.size() != B.superlattice.size()) {
    return superlattice.size() < B.superlattice.size();
  }
  return superlattice.superlattice() < B.superlattice.superlattice();
}

/// \brief Equality comparison of Supercell
bool Supercell::eq_impl(Supercell const &B) const {
  if (this == &B) {
    return true;
  }
  if (prim != B.prim) {
    throw std::runtime_error(
        "Error using Supercell::operator==(const Supercell& B): "
        "Only Supercell with shared prim may be compared this way.");
  }
  return superlattice.transformation_matrix_to_super() ==
         B.superlattice.transformation_matrix_to_super();
}

}  // namespace config
}  // namespace CASM
