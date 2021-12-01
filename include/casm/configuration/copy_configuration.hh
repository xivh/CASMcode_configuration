#ifndef CASM_config_copy_configuration
#define CASM_config_copy_configuration

#include "casm/configuration/definitions.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace config {

struct Configuration;
struct Supercell;

/// \brief Copy configuration DoF values into a supercell
Configuration copy_configuration(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell,
    UnitCell const &origin = UnitCell(0, 0, 0));

/// \brief Copy transformed configuration DoF values into a supercell
Configuration copy_configuration(
    Index factor_group_index, UnitCell translation, Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell,
    UnitCell const &origin = UnitCell(0, 0, 0));

/// \brief Return true if no translations within the supercell result in the
///     same configuration
bool is_primitive(Configuration const &configuration);

/// \brief Return the primitive configuration
Configuration make_primitive(Configuration const &configuration);

/// \brief Return the canonical configuration in the canonical supercell
Configuration make_in_canonical_supercell(Configuration const &configuration);

}  // namespace config
}  // namespace CASM

#endif
