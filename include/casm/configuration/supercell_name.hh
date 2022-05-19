#ifndef CASM_config_supercell_name
#define CASM_config_supercell_name

#include <string>
#include <vector>

namespace CASM {
namespace xtal {
class Lattice;
struct SymOp;
}  // namespace xtal
namespace config {

/// \brief Make the supercell name of a superlattice
std::string make_supercell_name(xtal::Lattice const &prim_lattice,
                                xtal::Lattice const &superlattice);

/// \brief Construct a superlattice from the supercell name
xtal::Lattice make_superlattice_from_supercell_name(
    xtal::Lattice const &prim_lattice, std::string supercell_name);

}  // namespace config
}  // namespace CASM

#endif
