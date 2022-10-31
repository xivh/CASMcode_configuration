#include "casm/configuration/ConfigurationSet.hh"

#include "casm/configuration/supercell_name.hh"

namespace CASM {
namespace config {

ConfigurationRecord::ConfigurationRecord(Configuration const &_configuration,
                                         std::string _supercell_name,
                                         std::string _configuration_id)
    : configuration(_configuration),
      supercell_name(_supercell_name),
      configuration_id(_configuration_id),
      configuration_name(supercell_name + "/" + configuration_id) {}

ConfigurationSet::ConfigurationSet(std::map<std::string, Index> _next_config_id)
    : m_next_config_id(_next_config_id) {}

bool ConfigurationSet::empty() const { return m_data.empty(); }

ConfigurationSet::size_type ConfigurationSet::size() const {
  return m_data.size();
}

void ConfigurationSet::clear() { return m_data.clear(); }

ConfigurationSet::const_iterator ConfigurationSet::begin() const {
  return m_data.begin();
}

ConfigurationSet::const_iterator ConfigurationSet::end() const {
  return m_data.end();
}

/// \brief Insert Configuration, setting supercell_name and
///     configuration_id automatically
std::pair<ConfigurationSet::iterator, bool> ConfigurationSet::insert(
    Configuration const &configuration) {
  auto const &superlattice = configuration.supercell->superlattice;
  std::string supercell_name = make_supercell_name(superlattice.prim_lattice(),
                                                   superlattice.superlattice());
  return this->insert(supercell_name, configuration);
}

/// \brief Insert Configuration with known supercell_name, setting
///     configuration_id automatically
std::pair<ConfigurationSet::iterator, bool> ConfigurationSet::insert(
    std::string const &supercell_name, Configuration const &configuration) {
  auto it = m_next_config_id.find(supercell_name);
  if (it == m_next_config_id.end()) {
    it = m_next_config_id.emplace(supercell_name, 0).first;
  }
  Index &configuration_id = it->second;

  auto res = m_data.insert(ConfigurationRecord(
      configuration, supercell_name, std::to_string(configuration_id)));
  if (res.second) {
    ++configuration_id;
  }
  return res;
}

/// \brief Insert ConfigurationRecord, allowing custom configuration_id
std::pair<ConfigurationSet::iterator, bool> ConfigurationSet::insert(
    ConfigurationRecord const &record) {
  return m_data.insert(record);
}

ConfigurationSet::const_iterator ConfigurationSet::find(
    Configuration const &configuration) const {
  ConfigurationRecord record(configuration, "", "");
  return m_data.find(record);
}

ConfigurationSet::const_iterator ConfigurationSet::find_by_name(
    std::string configuration_name) const {
  auto it = this->begin();
  auto end = this->end();
  for (; it != end; ++it) {
    if (it->configuration_name == configuration_name) {
      return it;
    }
  }
  return end;
}

ConfigurationSet::size_type ConfigurationSet::count(
    Configuration const &configuration) const {
  if (find(configuration) != end()) {
    return 1;
  }
  return 0;
}

ConfigurationSet::size_type ConfigurationSet::count_by_name(
    std::string configuration_name) const {
  if (find_by_name(configuration_name) != end()) {
    return 1;
  }
  return 0;
}

ConfigurationSet::const_iterator ConfigurationSet::erase(const_iterator it) {
  return m_data.erase(it);
}

ConfigurationSet::size_type ConfigurationSet::erase(
    Configuration const &configuration) {
  auto it = find(configuration);
  if (it == end()) {
    return 0;
  }
  m_data.erase(it);
  return 1;
}

ConfigurationSet::size_type ConfigurationSet::erase_by_name(
    std::string configuration_name) {
  auto it = find_by_name(configuration_name);
  if (it == end()) {
    return 0;
  }
  m_data.erase(it);
  return 1;
}

/// \brief Set IDs, by supercell_name, used to automatically ID new
/// configurations
void ConfigurationSet::set_next_config_id(
    std::map<std::string, Index> const &next_config_id) {
  m_next_config_id = next_config_id;
}

/// \brief IDs, by supercell_name, used to automatically ID new configurations
std::map<std::string, Index> const &ConfigurationSet::next_config_id() const {
  return m_next_config_id;
}

std::set<ConfigurationRecord> &ConfigurationSet::data() { return m_data; }

std::set<ConfigurationRecord> const &ConfigurationSet::data() const {
  return m_data;
}

/// \brief Make a map for finding ConfigurationRecord by configuration_name
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
}  // namespace CASM
