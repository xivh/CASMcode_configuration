#include "casm/configuration/io/json/Configuration_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clexulator/io/json/ConfigDoFValues_json_io.hh"
#include "casm/configuration/io/json/Supercell_json_io.hh"
#include "casm/configuration/supercell_name.hh"
#include "casm/misc/Validator.hh"

namespace CASM {

namespace {  // (anonymous)

/// \brief Validate ConfigDoFValues for correct DoF types and dimensions
///
/// \param validator An object to store error messages. May be a
///     standalone Validator, or an InputParser.
/// \param dof_values The ConfigDoFValues being validated
/// \param supercell_volume The integer supercell volume, as a multiple of the
/// prim \param basis_size The prim basis size \param global_dof_info The prim
/// global DoF basis sets \param local_dof_info the prim local DoF basis sets
///
/// Note:
/// - This does not validate that DoF values lie in the allowed DoF space
/// - This method could be moved to clexulator
void validate_dof_values(
    Validator &validator, clexulator::ConfigDoFValues const &dof_values,
    Index supercell_volume, Index basis_size,
    std::map<DoFKey, xtal::DoFSet> const &global_dof_info,
    std::map<DoFKey, std::vector<xtal::SiteDoFSet>> const &local_dof_info) {
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
    Index standard_dim =
        local_dof_info.at(name_value.first).front().basis().rows();
    if (name_value.second.rows() != standard_dim) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of local dof "
          << name_value.first
          << " is inconsistent with the standard dof dimension.";
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
    Index standard_dim = global_dof_info.at(name_value.first).basis().rows();
    if (name_value.second.rows() != standard_dim) {
      std::stringstream msg;
      msg << "Error reading Configuration from JSON: size of global dof "
          << name_value.first
          << " is inconsistent with the standard dof dimension.";
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

}  // namespace

namespace config {

ConfigurationRecord::ConfigurationRecord(Configuration const &_configuration,
                                         std::string _supercell_name,
                                         std::string _configuration_id,
                                         jsonParser const &_source,
                                         jsonParser const &_cache)
    : configuration(_configuration),
      supercell_name(_supercell_name),
      configuration_id(_configuration_id),
      source(_source),
      cache(_cache) {
  if (supercell_name.empty()) {
    supercell_name = make_supercell_name(
        configuration.supercell->superlattice.prim_lattice(),
        configuration.supercell->superlattice.superlattice());
  }
  configuration_name = supercell_name + "/" + configuration_id;
}

std::map<std::string, ConfigurationRecord const *>
make_index_by_configuration_name(
    std::set<ConfigurationRecord> const &configurations) {
  std::map<std::string, ConfigurationRecord const *> result;
  for (auto const &c : configurations) {
    result.emplace(c.configuration_name, &c);
  }
  return result;
}

}  // namespace config

void from_json(std::set<config::SupercellRecord> &supercells,
               std::set<config::ConfigurationRecord> &configurations,
               std::map<std::string, Index> &config_id, jsonParser const &json,
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

  if (!validator.valid()) {
    log.indent() << "Errors reading configurations:" << std::endl;
    report_and_throw_if_invalid(validator, log, error_if_invalid);
  }

  std::map<std::string, config::SupercellRecord const *>
      index_by_supercell_name = make_index_by_supercell_name(supercells);

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
      s = find_or_add_supercell_by_name(scel_it.name(), supercells,
                                        index_by_supercell_name, prim);
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
      from_json(dof_values, (*config_it)["dof"]);
      validate_dof_values(validator, dof_values,
                          s->supercell->unitcell_index_converter.total_sites(),
                          prim->basicstructure->basis().size(),
                          prim->global_dof_info, prim->local_dof_info);

      jsonParser source;
      config_it->get_if(source, "source");

      jsonParser cache;
      config_it->get_if(cache, "cache");

      if (!validator.valid()) {
        log.indent() << "Errors reading configurations:" << std::endl;
        log.indent() << "Error reading configuration: " << scel_it.name() << "/"
                     << config_it.name() << std::endl;
        report_and_throw_if_invalid(validator, log, error_if_invalid);
      }

      configurations.emplace(config::Configuration(s->supercell, dof_values),
                             scel_it.name(), config_it.name(), source, cache);
    }
  }

  // read next config id for each supercell
  config_id.clear();
  from_json(config_id, json["config_id"]);
}

jsonParser &to_json(std::set<config::ConfigurationRecord> const &configurations,
                    std::map<std::string, Index> const &config_id,
                    jsonParser &json) {
  json.put_obj();
  json["version"] = "1.0";
  json["supercells"] = jsonParser::object();
  for (const auto &c : configurations) {
    jsonParser &configjson =
        json["supercells"][c.supercell_name][c.configuration_id];
    to_json(c.configuration.dof_values, configjson["dof"]);

    if (c.source.size()) {
      to_json(c.source, configjson["source"]);
    }

    if (c.cache.size()) {
      to_json(c.cache, configjson["cache"]);
    }
  }

  json["config_id"] = config_id;
  return json;
}

config::Configuration jsonConstructor<config::Configuration>::from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  return std::move(
      *jsonMake<config::Configuration>::make_from_json(json, prim));
}

std::unique_ptr<config::Configuration>
jsonMake<config::Configuration>::make_from_json(
    jsonParser const &json, std::shared_ptr<config::Prim const> const &prim) {
  auto &log = CASM::log();
  ParentInputParser parser{json};
  std::runtime_error error_if_invalid{
      "Error reading Configuration from JSON input"};

  Eigen::Matrix3l T;
  parser.require(T, "transformation_matrix_to_supercell");
  report_and_throw_if_invalid(parser, log, error_if_invalid);
  auto supercell = std::make_shared<config::Supercell const>(prim, T);

  clexulator::ConfigDoFValues dof_values;
  parser.require(dof_values, "dof");

  validate_dof_values(parser, dof_values,
                      supercell->unitcell_index_converter.total_sites(),
                      prim->basicstructure->basis().size(),
                      prim->global_dof_info, prim->local_dof_info);

  report_and_throw_if_invalid(parser, log, error_if_invalid);

  return notstd::make_unique<config::Configuration>(supercell, dof_values);
}

/// Insert Configuration to JSON
jsonParser &to_json(config::Configuration const &configuration,
                    jsonParser &json);

/// Parser Configuration from JSON with error messages
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

  validate_dof_values(parser, dof_values,
                      supercell->unitcell_index_converter.total_sites(),
                      prim->basicstructure->basis().size(),
                      prim->global_dof_info, prim->local_dof_info);

  if (parser.valid()) {
    parser.value =
        notstd::make_unique<config::Configuration>(supercell, dof_values);
  }
}

/// Read Configuration from JSON
void from_json(config::Configuration &configuration, jsonParser const &json);

}  // namespace CASM
