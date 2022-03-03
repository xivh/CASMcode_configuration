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

/// \brief Make the canonical supercell name of a superlattice
std::string make_canonical_supercell_name(
    std::vector<xtal::SymOp> const &point_group,
    xtal::Lattice const &prim_lattice, xtal::Lattice const &superlattice);

/// \brief Construct a superlattice, in canonical form, from the
///     canonical supercell name
xtal::Lattice make_superlattice_from_canonical_supercell_name(
    std::vector<xtal::SymOp> const &point_group,
    xtal::Lattice const &prim_lattice, std::string canonical_supercell_name);

}  // namespace config
}  // namespace CASM

#endif
