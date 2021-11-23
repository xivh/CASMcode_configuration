#include "casm/configuration/Supercell.hh"

#include "casm/configuration/SupercellSymInfo.hh"

namespace CASM {
namespace config {

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Superlattice const &_superlattice)
    : prim(_prim),
      superlattice(_superlattice),
      unitcell_index_converter(superlattice.transformation_matrix_to_super()),
      unitcellcoord_index_converter(
          superlattice.transformation_matrix_to_super(),
          prim->basicstructure.basis().size()),
      sym_info(prim, superlattice, unitcell_index_converter,
               unitcellcoord_index_converter) {}

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Eigen::Matrix3l const &_superlattice_matrix)
    : Supercell(_prim, Superlattice(_prim->basicstructure.lattice(),
                                    _superlattice_matrix)) {}

}  // namespace config
}  // namespace CASM
