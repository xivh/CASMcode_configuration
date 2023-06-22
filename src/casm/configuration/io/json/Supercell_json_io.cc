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
      "Error reading Supercell from JSON input"};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  supercell = std::make_shared<config::Supercell const>(prim, T);
}

void from_json(std::shared_ptr<config::Supercell const> &supercell,
               jsonParser const &json, config::SupercellSet &supercells) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Supercell into SupercellSet from JSON input"};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  supercell = supercells.insert(T).first->supercell;
}

jsonParser &to_json(std::shared_ptr<config::Supercell const> const &supercell,
                    jsonParser &json) {
  if (!json.is_obj()) {
    throw std::runtime_error(
        "Error inserting supercell to json: not an object");
  }
  auto const &superlattice = supercell->superlattice;
  std::string supercell_name = config::make_supercell_name(
      superlattice.prim_lattice(), superlattice.superlattice());
  json["supercell_name"] = supercell_name;
  json["transformation_matrix_to_supercell"] =
      superlattice.transformation_matrix_to_super();
  return json;
}

/// \brief Read SupercellSet from JSON
void from_json(config::SupercellSet &supercells, jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim) {
  supercells.clear();

  std::set<std::string> matching_versions = {"1.0", "2.0"};

  // check json version
  if (!json.contains("version") ||
      !matching_versions.count(json["version"].get<std::string>())) {
    throw std::runtime_error(
        std::string("Error jsonDB incompatible version: found: ") +
        json["version"].get<std::string>());
  }

  if (!json.is_obj() || !json.contains("supercells")) {
    throw std::runtime_error("Error reading supercells: invalid format");
  }

  if (json.contains("supercells")) {
    auto it = json["supercells"].begin();
    auto end = json["supercells"].end();
    for (; it != end; ++it) {
      Eigen::Matrix3l mat;
      from_json(mat, *it);
      supercells.insert(std::make_shared<config::Supercell const>(prim, mat));
    }
  }
  if (json.contains("non_canonical_supercells")) {
    auto it = json["non_canonical_supercells"].begin();
    auto end = json["non_canonical_supercells"].end();
    for (; it != end; ++it) {
      if (!it->contains("transformation_matrix_to_supercell")) {
        throw std::runtime_error(
            "Error reading supercells: missing "
            "\"transformation_matrix_to_supercell\"");
      }
      Eigen::Matrix3l mat;
      from_json(mat, (*it)["transformation_matrix_to_supercell"]);
      supercells.insert(std::make_shared<config::Supercell const>(prim, mat));
    }
  }
}

/// \brief Write SupercellSet to JSON (version 2.0)
///
/// Notes:
/// - A version 1.0 file does not include non_canonical_supercells
///
/// Format:
/// \code
/// {
///   "version": "2.0",
///   "supercells": { // canonical supercells
///     <supercell_name>: <transformation_matrix_to_super>,
///     ...
///   },
///   "non_canonical_supercells": [
///     {
///       "supercell_name": <supercell_name>,
///       "canonical_supercell_name": <canonical_supercell_name>,
///       "transformation_matrix_to_supercell": <transformation_matrix_to_super>
///     },
///     ...
///   ]
/// }
/// \endcode
jsonParser &to_json(config::SupercellSet const &supercells, jsonParser &json) {
  json.put_obj();
  json["supercells"].put_obj();
  json["non_canonical_supercells"].put_array();
  json["version"] = "2.0";
  for (auto const &s : supercells) {
    if (s.is_canonical) {
      json["supercells"][s.supercell_name] =
          s.supercell->superlattice.transformation_matrix_to_super();
    } else {
      jsonParser tjson;
      tjson["canonical_supercell_name"] = s.canonical_supercell_name;
      tjson["supercell_name"] = s.supercell_name;
      tjson["transformation_matrix_to_supercell"] =
          s.supercell->superlattice.transformation_matrix_to_super();
      json["non_canonical_supercells"].push_back(tjson);
    }
  }
  return json;
}

/// \brief Write SupercellSet to JSON (version 1.0)
///
/// Format:
/// \code
/// {
///   "version": "1.0",
///   "supercells": { // canonical supercells
///     <supercell_name>: <transformation_matrix_to_super>,
///     ...
///   }
/// }
/// \endcode
jsonParser &to_json_v1(config::SupercellSet const &supercells,
                       jsonParser &json) {
  json.put_obj();
  json["supercells"].put_obj();
  json["version"] = "1.0";
  for (auto const &s : supercells) {
    if (s.is_canonical) {
      json["supercells"][s.supercell_name] =
          s.supercell->superlattice.transformation_matrix_to_super();
    }
  }
  return json;
}

}  // namespace CASM
