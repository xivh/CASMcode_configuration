#ifndef CASM_config_Supercell_json_io
#define CASM_config_Supercell_json_io

#include <memory>
#include <set>

namespace CASM {

class jsonParser;

namespace config {
struct Prim;
struct Supercell;
struct SupercellRecord;
class SupercellSet;
}  // namespace config

void from_json(std::shared_ptr<config::Supercell const> &supercell,
               jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim);

void from_json(std::shared_ptr<config::Supercell const> &supercell,
               jsonParser const &json, config::SupercellSet &supercells);

jsonParser &to_json(std::shared_ptr<config::Supercell const> const &supercell,
                    jsonParser &json);

/// \brief Read SupercellSet from JSON
void from_json(config::SupercellSet &supercells, jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim);

/// \brief Write SupercellSet to JSON (version 2.0)
jsonParser &to_json(config::SupercellSet const &supercells, jsonParser &json);

/// \brief Write SupercellSet to JSON (version 1.0)
jsonParser &to_json_v1(config::SupercellSet const &supercells,
                       jsonParser &json);

}  // namespace CASM

#endif
