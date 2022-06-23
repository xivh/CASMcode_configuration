#ifndef CASM_config_Configuration_json_io
#define CASM_config_Configuration_json_io

#include <map>
#include <set>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;

namespace config {

struct SupercellRecord;

struct ConfigurationRecord {
  ConfigurationRecord(Configuration const &_configuration,
                      std::string _supercell_name = std::string(),
                      std::string _configuration_id = std::string("none"),
                      jsonParser const &_source = jsonParser(),
                      jsonParser const &_cache = jsonParser());

  /// \brief Shared pointer to the configuration
  Configuration configuration;

  /// \brief Name of supercell for the configuration (i.e. "SCEL4_2_2_1_0_0_0")
  std::string supercell_name;

  /// \brief Distinguish configurations in the same supercell (i.e. "2")
  std::string configuration_id;

  /// \brief Supercell name and configuration id (i.e. "SCEL4_2_2_1_0_0_0/2")
  std::string configuration_name;

  /// \brief Optional, can hold general information about how the configuration
  /// was generated
  jsonParser source;

  /// \brief Optional, can hold cached data about the configuration that is
  /// expensive to calculate
  jsonParser cache;

  bool operator<(ConfigurationRecord const &rhs) const {
    return this->configuration < rhs.configuration;
  }
};

std::map<std::string, ConfigurationRecord const *>
make_index_by_configuration_name(
    std::set<ConfigurationRecord> const &configurations);

}  // namespace config

void from_json(std::set<config::SupercellRecord> &supercells,
               std::set<config::ConfigurationRecord> &configurations,
               std::map<std::string, Index> &config_id, jsonParser const &json,
               std::shared_ptr<config::Prim const> const &prim);

jsonParser &to_json(std::set<config::ConfigurationRecord> const &configurations,
                    std::map<std::string, Index> const &config_id,
                    jsonParser &json);

template <typename T>
struct jsonConstructor;
template <typename T>
struct jsonMake;
class jsonParser;

template <>
struct jsonConstructor<config::Configuration> {
  /// Read Configuration from JSON
  static config::Configuration from_json(
      jsonParser const &json, std::shared_ptr<config::Prim const> const &prim);
};

template <>
struct jsonMake<config::Configuration> {
  /// Read Configuration from JSON
  static std::unique_ptr<config::Configuration> make_from_json(
      jsonParser const &json, std::shared_ptr<config::Prim const> const &prim);
};

/// Insert Configuration to JSON
jsonParser &to_json(config::Configuration const &configuration,
                    jsonParser &json);

/// Parser Configuration from JSON with error messages
void parse(InputParser<config::Configuration> &parser,
           std::shared_ptr<config::Prim const> const &prim);

/// Read Configuration from JSON
void from_json(config::Configuration &configuration, jsonParser const &json);

}  // namespace CASM

#endif
