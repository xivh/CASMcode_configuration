#include "casm/configuration/io/json/Supercell_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/supercell_name.hh"

namespace CASM {

void from_json(std::shared_ptr<config::Supercell const> &supercell,
               jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Configuration from JSON input"};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  supercell = std::make_shared<config::Supercell const>(prim, T);
}

jsonParser &to_json(std::shared_ptr<config::Supercell const> const &supercell,
                    jsonParser &json);

void from_json(config::SupercellSet &supercells, jsonParser const &json,
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
    supercells.insert(std::make_shared<config::Supercell const>(prim, mat));
  }
}

jsonParser &to_json(config::SupercellSet const &supercells, jsonParser &json) {
  json.put_obj();
  json["version"] = "1.0";
  for (auto const &s : supercells) {
    json["supercells"][s.supercell_name] =
        s.supercell->superlattice.transformation_matrix_to_super();
  }
  return json;
}

}  // namespace CASM
