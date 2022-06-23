#ifndef CASM_config_Supercell_json_io
#define CASM_config_Supercell_json_io

#include <map>
#include <set>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure for holding / reading / writing supercells
struct SupercellRecord {
  SupercellRecord(std::shared_ptr<Supercell const> const &_supercell);

  std::shared_ptr<config::Supercell const> supercell;

  std::string supercell_name;

  bool operator<(SupercellRecord const &rhs) const {
    return *this->supercell < *rhs.supercell;
  }
};

/// \brief Make a map for finding SupercellRecord by supercell_name
std::map<std::string, SupercellRecord const *> make_index_by_supercell_name(
    std::set<SupercellRecord> const &supercells);

/// \brief Find or add supercell by name
SupercellRecord const *find_or_add_supercell_by_name(
    std::string const &supercell_name,
    std::set<config::SupercellRecord> &supercells,
    std::map<std::string, SupercellRecord const *> &index_by_supercell_name,
    std::shared_ptr<config::Prim const> const &prim);

}  // namespace config

void from_json(std::set<config::SupercellRecord> &supercells,
               jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim);

jsonParser &to_json(std::set<config::SupercellRecord> const &supercells,
                    jsonParser &json);

}  // namespace CASM

#endif
