#include "casm/configuration/io/json/Supercell_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/configuration/supercell_name.hh"

namespace CASM {
namespace config {

SupercellRecord::SupercellRecord(
    std::shared_ptr<Supercell const> const &_supercell)
    : supercell(throw_if_equal_to_nullptr(
          _supercell,
          "Error in SupercellRecord constructor: value == nullptr")),
      supercell_name(
          make_supercell_name(supercell->superlattice.prim_lattice(),
                              supercell->superlattice.superlattice())) {}

std::map<std::string, SupercellRecord const *> make_index_by_supercell_name(
    std::set<SupercellRecord> const &supercells) {
  std::map<std::string, SupercellRecord const *> result;
  for (auto const &s : supercells) {
    result.emplace(s.supercell_name, &s);
  }
  return result;
}

/// \brief Find or add supercell by name
SupercellRecord const *find_or_add_supercell_by_name(
    std::string const &supercell_name,
    std::set<config::SupercellRecord> &supercells,
    std::map<std::string, SupercellRecord const *> &index_by_supercell_name,
    std::shared_ptr<config::Prim const> const &prim) {
  SupercellRecord const *s = nullptr;
  auto it = index_by_supercell_name.find(supercell_name);
  if (it == index_by_supercell_name.end()) {
    auto supercell = std::make_shared<Supercell const>(
        prim, make_superlattice_from_supercell_name(
                  prim->basicstructure->lattice(), supercell_name));
    s = &*supercells.insert(supercell).first;
    if (supercell_name != s->supercell_name) {
      throw std::runtime_error(
          "Error in find_or_add_supercell_by_name: supercell_name mismatch");
    }
    index_by_supercell_name.emplace(supercell_name, s);
  } else {
    s = it->second;
  }
  return s;
}

}  // namespace config

void from_json(std::set<config::SupercellRecord> &supercells,
               jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim) {
  supercells.clear();

  // check json version
  if (!json.contains("version") ||
      json["version"].get<std::string>() != "1.0") {
    throw std::runtime_error(
        std::string("Error jsonDB version mismatch: found: ") +
        json["version"].get<std::string>() + " expected: 1.0");
  }

  if (!json.is_obj() || !json.contains("supercells")) {
    throw std::runtime_error("Error reading supercells: invalid format");
  }

  auto it = json["supercells"].begin();
  auto end = json["supercells"].end();
  for (; it != end; ++it) {
    Eigen::Matrix3l mat;
    from_json(mat, *it);
    supercells.emplace(std::make_shared<config::Supercell const>(prim, mat));
  }
}

jsonParser &to_json(std::set<config::SupercellRecord> const &supercells,
                    jsonParser &json) {
  json.put_obj();
  json["version"] = "1.0";
  for (auto const &s : supercells) {
    json["supercells"][s.supercell_name] =
        s.supercell->superlattice.transformation_matrix_to_super();
  }
  return json;
}

}  // namespace CASM
