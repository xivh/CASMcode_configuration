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

}  // namespace config
}  // namespace CASM
