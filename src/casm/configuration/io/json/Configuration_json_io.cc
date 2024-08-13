#include "casm/configuration/io/json/Configuration_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/configuration/ConfigurationSet.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/io/json/Supercell_json_io.hh"
#include "casm/configuration/supercell_name.hh"
#include "casm/misc/Validator.hh"

namespace CASM {

namespace {  // (anonymous)

/// \brief Validate ConfigDoFValues for correct DoF types and dimensions for
///     the standard basis
///
/// \param validator An object to store error messages. May be a
///     standalone Validator, or an InputParser.
/// \param dof_values The ConfigDoFValues being validated
/// \param supercell_volume The integer supercell volume, as a multiple of the
///     prim
/// \param basis_size The prim basis size
/// \param global_dof_info The prim global DoF basis sets
/// \param local_dof_info the prim local DoF basis sets
/// \param standard_basis If true, expect DoF values to have the dimension of
/// the standard basis. Otherwise, expect DoF values to have the dimension of
/// the prim basis.
///
/// Note:
/// - This does not validate that DoF values lie in the allowed DoF space
/// - This method could be moved to clexulator
void validate_dof_values(
    Validator &validator, clexulator::ConfigDoFValues const &dof_values,
    Index supercell_volume, Index basis_size,
    std::map<DoFKey, xtal::DoFSet> const &global_dof_info,
    std::map<DoFKey, std::vector<xtal::SiteDoFSet>> const &local_dof_info,
    bool standard_basis) {
  Index N_site = supercell_volume * basis_size;

  if (dof_values.occupation.size() != N_site) {
    std::stringstream msg;
    msg << "Error reading Configuration from JSON: "
        << "supercell number of sites (" << N_site << ") and "
        << "occupation size (" << dof_values.occupation.size()
        << ") are inconsistent.";
    validator.error.insert(msg.str());
  }

  // check for missing local dof
  for (auto const &name_info : local_dof_info) {
    if (!dof_values.local_dof_values.count(name_info.first)) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: local dof "
          << name_info.first << " is missing.";
      validator.error.insert(msg.str());
    }
  }

  // check local dof type and dimension
  for (auto const &name_value : dof_values.local_dof_values) {
    if (!local_dof_info.count(name_value.first)) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: local dof "
          << name_value.first << " is inconsistent with the prim.";
      validator.error.insert(msg.str());
      continue;
    }
    if (name_value.second.cols() != N_site) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of local dof "
          << name_value.first
          << " is inconsistent with the supercell number of sites.";
      validator.error.insert(msg.str());
    }
    Index expected_dim;
    if (standard_basis) {
      expected_dim = local_dof_info.at(name_value.first).front().basis().rows();
    } else {
      expected_dim = clexulator::max_dim(local_dof_info.at(name_value.first));
    }
    if (name_value.second.rows() != expected_dim) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of local dof "
          << name_value.first << " is inconsistent with the "
          << (standard_basis ? "standard basis" : "max prim basis")
          << " dof dimension.";
      validator.error.insert(msg.str());
    }
  }

  // check for missing global dof
  for (auto const &name_info : global_dof_info) {
    if (!dof_values.global_dof_values.count(name_info.first)) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: global dof "
          << name_info.first << " is missing.";
      validator.error.insert(msg.str());
    }
  }

  // check global dof type and dimension
  for (auto const &name_value : dof_values.global_dof_values) {
    if (!global_dof_info.count(name_value.first)) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: global dof "
          << name_value.first << " is inconsistent with the prim.";
      validator.error.insert(msg.str());
      continue;
    }
    if (name_value.second.cols() != 1) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of global dof "
          << name_value.first << " is incorrect.";
      validator.error.insert(msg.str());
    }
    Index expected_dim;
    if (standard_basis) {
      expected_dim = global_dof_info.at(name_value.first).basis().rows();
    } else {
      expected_dim = global_dof_info.at(name_value.first).basis().cols();
    }
    if (name_value.second.rows() != expected_dim) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of global dof "
          << name_value.first << " is inconsistent with the "
          << (standard_basis ? "standard basis" : "max prim basis")
          << " dof dimension.";
      validator.error.insert(msg.str());
    }
  }
}

template <typename ErrorType>
void report_and_throw_if_invalid(Validator const &validator, Log &log,
                                 ErrorType error_if_invalid) {
  if (!validator.valid()) {
    log.increase_indent();
    for (auto const &e : validator.error) {
      log.indent() << e << std::endl;
    }
    log.decrease_indent();
    throw error_if_invalid;
  }
}

std::map<std::string, Eigen::MatrixXd> _empty_local() {
  return std::map<std::string, Eigen::MatrixXd>();
}

std::map<std::string, Eigen::VectorXd> _empty_global() {
  return std::map<std::string, Eigen::VectorXd>();
}

std::unique_ptr<config::ConfigurationWithProperties>
_make_unique_with_properties(config::Configuration const &config) {
  return notstd::make_unique<config::ConfigurationWithProperties>(
      config, _empty_local(), _empty_global());
}

/// Parser properties for ConfigurationWithProperties from JSON with error
/// messages
///
/// Requires parser.value != nullptr
void _parse_properties(
    InputParser<config::ConfigurationWithProperties> &parser) {
  if (parser.value == nullptr) {
    throw std::runtime_error(
        "Error parsing properties for ConfigurationWithProperties");
  }
  config::ConfigurationWithProperties &self = *parser.value;

  // local properties
  auto json_it = parser.self.find("local_properties");
  if (json_it != parser.self.end()) {
    for (auto it = json_it->begin(); it != json_it->end(); ++it) {
      Eigen::MatrixXd tvalues;
      fs::path values_path =
          fs::path{"local_properties"} / it.name() / "values";
      parser.require(tvalues, values_path);
      try {
        self.local_properties.emplace(it.name(), tvalues.transpose());
      } catch (std::exception const &e) {
        parser.insert_error(values_path, e.what());
      }
    }
  }

  // global properties
  json_it = parser.self.find("global_properties");
  if (json_it != parser.self.end()) {
    for (auto it = json_it->begin(); it != json_it->end(); ++it) {
      Eigen::VectorXd tvalues;
      fs::path values_path =
          fs::path{"global_properties"} / it.name() / "values";
      parser.require(tvalues, values_path);
      try {
        self.global_properties.emplace(it.name(), tvalues.transpose());
      } catch (std::exception const &e) {
        parser.insert_error(values_path, e.what());
      }
    }
  }
}

/// \brief Check if JSON-formatted Configuration DoF values are in the prim
/// basis
///
/// Notes:
/// - If "basis": "prim", DoF values are in the prim basis -> return true
/// - If "basis": "standard", DoF values are in the standard basis -> return
/// false
/// - If "basis" is something else, return false and insert an error message
/// - If "basis" is not present, assume DoF values are in the standard basis ->
/// return false
///
/// \returns True, if the Configuration DoF values stored in `json` are
///     expressed in the prim basis; return false otherwwise.
bool get_read_prim_basis(Validator &validator, jsonParser const &json,
                         std::string what) {
  bool read_prim_basis = false;
  if (json.contains("basis")) {
    if (json["basis"].is_string() &&
        json["basis"].get<std::string>() == "prim") {
      read_prim_basis = true;
    } else if (json["basis"].is_string() &&
               json["basis"].get<std::string>() == "standard") {
      read_prim_basis = false;
    } else {
      validator.error.insert("Error reading " + what +
                             ": If present, \"basis\" value must be "
                             "\"prim\" or \"standard\".");
    }
  }
  return read_prim_basis;
}

clexulator::ConfigDoFValues read_dof_values(
    Validator &validator, jsonParser const &json,
    std::shared_ptr<config::Supercell const> supercell, bool read_prim_basis) {
  clexulator::ConfigDoFValues dof_values;
  if (!json.contains("dof")) {
    validator.error.insert("Error reading DoF values: \"dof\" not found.");
    return dof_values;
  }
  auto &prim = *supercell->prim;
  from_json(dof_values, json["dof"]);
  validate_dof_values(validator, dof_values,
                      supercell->unitcell_index_converter.total_sites(),
                      prim.basicstructure->basis().size(), prim.global_dof_info,
                      prim.local_dof_info, !read_prim_basis);
  if (!read_prim_basis) {
    dof_values = clexulator::from_standard_values(
        dof_values, prim.basicstructure->basis().size(),
        supercell->unitcell_index_converter.total_sites(), prim.global_dof_info,
        prim.local_dof_info);
  }
  return dof_values;
}
}  // namespace

void from_json(config::SupercellSet &supercells,
               config::ConfigurationSet &configurations, jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim) {
  configurations.clear();

  auto &log = CASM::log();
  Validator validator;
  std::runtime_error error_if_invalid{"Error reading configurations"};

  // check json version
  if (!json.contains("version") ||
      json["version"].get<std::string>() != "1.0") {
    validator.error.insert(
        std::string("Error jsonDB version mismatch: found: ") +
        json["version"].get<std::string>() + " expected: 1.0");
  }

  if (!json.is_obj() || !json.contains("supercells")) {
    validator.error.insert("Error reading configurations: invalid format");
  }

  bool read_prim_basis =
      get_read_prim_basis(validator, json, "ConfigurationSet");

  if (!validator.valid()) {
    log.indent() << "Errors reading configurations:" << std::endl;
    report_and_throw_if_invalid(validator, log, error_if_invalid);
  }

  std::map<std::string, config::SupercellRecord const *>
      index_by_supercell_name =
          make_index_by_canonical_supercell_name(supercells.data());

  // read config list contents
  auto scel_it = json["supercells"].begin();
  auto scel_end = json["supercells"].end();

  clexulator::ConfigDoFValues dof_values;

  for (; scel_it != scel_end; ++scel_it) {
    auto config_it = scel_it->begin();
    auto config_end = scel_it->end();

    // try to find or add supercell by name
    config::SupercellRecord const *s = nullptr;
    try {
      s = find_or_add_canonical_supercell_by_name(
          scel_it.name(), supercells.data(), index_by_supercell_name, prim);
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error: could not find or construct supercell '" << scel_it.name()
          << "' by name: " << e.what();
      validator.error.insert(msg.str());
    }

    if (!validator.valid()) {
      log.indent() << "Error reading configurations:" << std::endl;
      report_and_throw_if_invalid(validator, log, error_if_invalid);
    }

    // try to construct configurations for supercell
    for (; config_it != config_end; ++config_it) {
      dof_values = read_dof_values(validator, (*config_it), s->supercell,
                                   read_prim_basis);

      // jsonParser source;
      // config_it->get_if(source, "source");
      //
      // jsonParser cache;
      // config_it->get_if(cache, "cache");

      if (!validator.valid()) {
        log.indent() << "Errors reading configurations:" << std::endl;
        log.indent() << "Error reading configuration: " << scel_it.name() << "/"
                     << config_it.name() << std::endl;
        report_and_throw_if_invalid(validator, log, error_if_invalid);
      }

      config::ConfigurationRecord record(
          config::Configuration(s->supercell, dof_values), scel_it.name(),
          config_it.name());
      configurations.insert(record);
    }
  }

  // read next config id for each supercell
  std::map<std::string, Index> next_config_id;
  from_json(next_config_id, json["config_id"]);
  configurations.set_next_config_id(next_config_id);
}

/// \brief Write ConfigurationSet to JSON
///
/// \param configurations The ConfigurationSet
/// \param json The jsonParser
/// \param write_prim_basis If true, write DoF values using the prim basis.
///     Default (false) is to write DoF values in the standard basis.
/// \return jsonParser with ConfigurationSet added.
jsonParser &to_json(config::ConfigurationSet const &configurations,
                    jsonParser &json, bool write_prim_basis) {
  json.put_obj();
  json["version"] = "1.0";
  json["supercells"] = jsonParser::object();
  if (write_prim_basis) {
    json["basis"] = "prim";
  } else {
    json["basis"] = "standard";
  }
  for (const auto &c : configurations) {
    jsonParser &configjson =
        json["supercells"][c.supercell_name][c.configuration_id];
    if (write_prim_basis) {
      to_json(c.configuration.dof_values, configjson["dof"]);
    } else {
      to_json(make_standard_dof_values(c.configuration), configjson["dof"]);
    }

    // if (c.source.size()) {
    //   to_json(c.source, configjson["source"]);
    // }
    //
    // if (c.cache.size()) {
    //   to_json(c.cache, configjson["cache"]);
    // }
  }

  json["config_id"] = configurations.next_config_id();
  return json;
}

config::Configuration jsonConstructor<config::Configuration>::from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  return std::move(
      *jsonMake<config::Configuration>::make_from_json(json, prim));
}

config::Configuration jsonConstructor<config::Configuration>::from_json(
    jsonParser const &json, config::SupercellSet &supercells) {
  return std::move(
      *jsonMake<config::Configuration>::make_from_json(json, supercells));
}

std::unique_ptr<config::Configuration>
jsonMake<config::Configuration>::make_from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Configuration from JSON input"};
  auto result = parser.parse_as<config::Configuration>(prim);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  return std::move(result->value);
}

std::unique_ptr<config::Configuration>
jsonMake<config::Configuration>::make_from_json(
    jsonParser const &json, config::SupercellSet &supercells) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Configuration from JSON input"};
  auto result = parser.parse_as<config::Configuration>(supercells);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  return std::move(result->value);
}

/// Insert Configuration to JSON
///
/// \param configuration A Configuration
/// \param json A jsonParser
/// \param write_prim_basis If true, write DoF values using the prim basis.
///     Default (false) is to write DoF values in the standard basis.
///     DoF values in `configuration` are expected to be in the prim basis.
///
jsonParser &to_json(config::Configuration const &configuration,
                    jsonParser &json, bool write_prim_basis) {
  if (!json.is_obj()) {
    throw std::runtime_error(
        "Error inserting configuration to json: not an object");
  }
  auto const &superlattice = configuration.supercell->superlattice;
  std::string supercell_name = config::make_supercell_name(
      superlattice.prim_lattice(), superlattice.superlattice());
  json["supercell_name"] = supercell_name;
  json["transformation_matrix_to_supercell"] =
      superlattice.transformation_matrix_to_super();
  if (write_prim_basis) {
    json["basis"] = "prim";
    to_json(configuration.dof_values, json["dof"]);
  } else {
    json["basis"] = "standard";
    to_json(make_standard_dof_values(configuration), json["dof"]);
  }

  return json;
}

/// Parse Configuration from JSON with error messages
///
/// Notes:
/// - This version creates a new Supercell for each Configuration.
///   Comparisons will still work correctly, but if reading many
///   Configuration, memory usage will be increased and speed will
///   be decreased.
void parse(InputParser<config::Configuration> &parser,
           std::shared_ptr<config::Prim const> const &prim) {
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  auto supercell = std::make_shared<config::Supercell const>(prim, T);

  clexulator::ConfigDoFValues dof_values;
  parser.require(dof_values, "dof");

  bool read_prim_basis =
      get_read_prim_basis(parser, parser.self, "Configuration");
  dof_values = read_dof_values(parser, parser.self, supercell, read_prim_basis);

  if (parser.valid()) {
    parser.value =
        notstd::make_unique<config::Configuration>(supercell, dof_values);
  }
}

/// Parse Configuration from JSON with error messages
///
/// Notes:
/// - This version avoids duplicate Supercells.
void parse(InputParser<config::Configuration> &parser,
           config::SupercellSet &supercells) {
  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  auto supercell = supercells.insert(T).first->supercell;

  clexulator::ConfigDoFValues dof_values;
  parser.require(dof_values, "dof");

  bool read_prim_basis =
      get_read_prim_basis(parser, parser.self, "Configuration");
  dof_values = read_dof_values(parser, parser.self, supercell, read_prim_basis);

  if (parser.valid()) {
    parser.value =
        notstd::make_unique<config::Configuration>(supercell, dof_values);
  }
}

// --- ConfigurationWithProperties

config::ConfigurationWithProperties
jsonConstructor<config::ConfigurationWithProperties>::from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  return std::move(
      *jsonMake<config::ConfigurationWithProperties>::make_from_json(json,
                                                                     prim));
}

config::ConfigurationWithProperties
jsonConstructor<config::ConfigurationWithProperties>::from_json(
    jsonParser const &json, config::SupercellSet &supercells) {
  return std::move(
      *jsonMake<config::ConfigurationWithProperties>::make_from_json(
          json, supercells));
}

std::unique_ptr<config::ConfigurationWithProperties>
jsonMake<config::ConfigurationWithProperties>::make_from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading ConfigurationWithProperties from JSON input"};
  auto result = parser.parse_as<config::ConfigurationWithProperties>(prim);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  return std::move(result->value);
}

std::unique_ptr<config::ConfigurationWithProperties>
jsonMake<config::ConfigurationWithProperties>::make_from_json(
    jsonParser const &json, config::SupercellSet &supercells) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading ConfigurationWithProperties from JSON input"};
  auto result =
      parser.parse_as<config::ConfigurationWithProperties>(supercells);
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  return std::move(result->value);
}

/// Insert ConfigurationWithProperties to JSON
///
/// \param configuration A Configuration
/// \param json A jsonParser
/// \param write_prim_basis If true, write DoF values using the prim basis.
///     Default (false) is to write DoF values in the standard basis.
///
/// Format:
/// \code
/// {
///   "configuration": <configuration object>,
///   "global_properties": {
///     <property_name>: {
///       "values": <1d array of numbers, property values as provided>
///     },
///     ...
///   },
///   "local_properties": {
///     <property_name>: {
///       "values": <2d array of numbers, property values as provided>
///     },
///     ...
///   }
/// }
/// \endcode
///
/// Reminder about standard DoF basis vs prim DoF basis:
/// - Properties are expected to be in the standard basis
/// - This function does not convert ConfigDoFValues between bases, it writes
///   values as they are.
/// - Conversions, if necessary, must be done before calling `to_json` / after
///   calling `from_json`.
///
/// Note:
/// - "configuration":
///       The configuration, as JSON
///
/// - "local_properties":
///       The local property values are represented as a matrix, with each row
///       representing a site property value and each column representing a
///       component of the property value:
///           number of cols = property dimension (i.e. 3 for "disp")
///           number of rows = (Supercell volume as multiple of the prim) *
///               (prim basis size).
///
///       Example: "disp" values, supercell volume=3, prim basis size=2
///
///           "local_properties": {
///             "disp": {
///               "values": [
///                 [dx[1], dy[1], dz[1]], // "disp" values on site: sublattice
///                 0, unit cell 0 [dx[2], dy[2], dz[2]], // "disp" values on
///                 site: sublattice 0, unit cell 1 [dx[3], dy[3], dz[3]], //
///                 "disp" values on site: sublattice 0, unit cell 2 [dx[4],
///                 dy[4], dz[4]], // "disp" values on site: sublattice 1, unit
///                 cell 0 [dx[5], dy[5], dz[5]], // "disp" values on site:
///                 sublattice 1, unit cell 1 [dx[6], dy[6], dz[6]], // "disp"
///                 values on site: sublattice 1, unit cell 2
///               ]
///             }
///           }
///
/// - "global_properties":
///       The global property values are represented as a vector of size equal
///       to the dimension of the DoF (i.e. 6 for "GLstrain" in standard basis).
///
///       Example: "GLstrain" values (any supercell volume and prim basis size)
///
///           "global_properties": {
///             "GLstrain": {
///               "values": [Exx, Eyy, Ezz, sqrt(2)Eyz, sqrt(2)Exz, sqrt(2)Exy]
///             }
///           }
///
jsonParser &to_json(
    config::ConfigurationWithProperties const &configuration_with_properties,
    jsonParser &json, bool write_prim_basis) {
  auto const &x = configuration_with_properties;
  if (!json.is_obj()) {
    throw std::runtime_error(
        "Error inserting ConfigurationWithProperties to json: not an object");
  }
  to_json(x.configuration, json["configuration"], write_prim_basis);
  if (!x.local_properties.empty()) {
    for (auto const &local_property : x.local_properties) {
      to_json(local_property.second.transpose(),
              json["local_properties"][local_property.first]["values"]);
    }
  }
  if (!x.global_properties.empty()) {
    for (auto const &global_property : x.global_properties) {
      to_json_array(global_property.second,
                    json["global_properties"][global_property.first]["values"]);
    }
  }
  return json;
}

/// Parse ConfigurationWithProperties from JSON with error messages
///
/// Notes:
/// - This version creates a new Supercell for each Configuration.
///   Comparisons will still work correctly, but if reading many
///   Configuration, memory usage will be increased and speed will
///   be decreased.
void parse(InputParser<config::ConfigurationWithProperties> &parser,
           std::shared_ptr<config::Prim const> const &prim) {
  // parse configuration
  std::unique_ptr<config::Configuration> config =
      parser.require<config::Configuration>("configuration", prim);
  if (config == nullptr) {
    return;
  }

  // parse local_properties and global_properties
  parser.value = _make_unique_with_properties(*config);
  _parse_properties(parser);
  if (!parser.valid()) {
    parser.value.reset();
  }
}

/// Parser ConfigurationWithProperties from JSON with error messages
///
/// Notes:
/// - This version avoids duplicate Supercells.
void parse(InputParser<config::ConfigurationWithProperties> &parser,
           config::SupercellSet &supercells) {
  // parse configuration
  std::unique_ptr<config::Configuration> config =
      parser.require<config::Configuration>("configuration", supercells);
  if (config == nullptr) {
    return;
  }

  // parse local_properties and global_properties
  parser.value = _make_unique_with_properties(*config);
  _parse_properties(parser);
  if (!parser.valid()) {
    parser.value.reset();
  }
}

}  // namespace CASM
