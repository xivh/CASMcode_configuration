#include "casm/configuration/Supercell.hh"

#include "casm/configuration/SupercellSymInfo.hh"

namespace CASM {
namespace config {

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Superlattice const &_superlattice)
    : prim(_prim),
      superlattice(_superlattice),
      unitcell_index_converter(_superlattice.transformation_matrix_to_super()),
      unitcellcoord_index_converter(
          _superlattice.transformation_matrix_to_super(),
          _prim->basicstructure.basis().size()),
      sym_info(_prim, _superlattice) {}

}  // namespace config
}  // namespace CASM
