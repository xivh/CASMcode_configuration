#ifndef CASM_config_ConfigurationSet
#define CASM_config_ConfigurationSet

#include <map>
#include <set>

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace config {

/// \brief Data structure for holding / reading / writing configurations
struct ConfigurationRecord : public Comparisons<CRTPBase<ConfigurationRecord>> {
  ConfigurationRecord(Configuration const &_configuration,
                      std::string _supercell_name,
                      std::string _configuration_id);

  /// \brief Shared pointer to the configuration
  Configuration configuration;

  /// \brief Name of canonical supercell for the configuration (i.e.
  /// "SCEL4_2_2_1_0_0_0")
  std::string supercell_name;

  /// \brief Distinguish configurations in the same canonical supercell (i.e.
  /// "2")
  std::string configuration_id;

  /// \brief Canonical supercell name and configuration id (i.e.
  /// "SCEL4_2_2_1_0_0_0/2")
  std::string configuration_name;

  bool operator<(ConfigurationRecord const &rhs) const {
    return this->configuration < rhs.configuration;
  }

 private:
  friend struct Comparisons<CRTPBase<ConfigurationRecord>>;
};

/// \brief Data structure for holding / reading / writing canonical
/// configurations
///
/// Notes:
/// - This class is **only** valid for use holding configurations with
///   canonical supercells
/// - It is valid for use with non-canonical configurations or non-primitive
///   configurations, if they are in the canonical supercell
/// - For configurations with non-canonical supercells, do not use this class,
///   use a std::vector<Configuration> or other container
/// - Includes a map of supercell_name -> next configuration id that can
///   be used to automatically provide new configurations with sequential IDs
class ConfigurationSet {
 public:
  ConfigurationSet(std::map<std::string, Index> _next_config_id = {});

  typedef std::set<ConfigurationRecord>::size_type size_type;
  typedef std::set<ConfigurationRecord>::iterator iterator;
  typedef std::set<ConfigurationRecord>::const_iterator const_iterator;

  bool empty() const;

  size_type size() const;

  void clear();

  const_iterator begin() const;

  const_iterator end() const;

  /// \brief Insert Configuration, setting supercell_name and
  ///     configuration_id automatically
  std::pair<iterator, bool> insert(Configuration const &configuration);

  /// \brief Insert Configuration with known supercell_name, setting
  ///     configuration_id automatically
  std::pair<iterator, bool> insert(std::string const &supercell_name,
                                   Configuration const &configuration);

  /// \brief Insert ConfigurationRecord, allowing custom configuration_id
  std::pair<iterator, bool> insert(ConfigurationRecord const &record);

  const_iterator find(Configuration const &configuration) const;

  const_iterator find_by_name(std::string configuration_name) const;

  size_type count(Configuration const &configuration) const;

  size_type count_by_name(std::string configuration_name) const;

  const_iterator erase(const_iterator it);

  size_type erase(Configuration const &configuration);

  size_type erase_by_name(std::string configuration_name);

  /// \brief Set IDs, by supercell_name, used to automatically ID new
  /// configurations
  void set_next_config_id(std::map<std::string, Index> const &next_config_id);

  /// \brief IDs, by supercell_name, used to automatically ID new configurations
  std::map<std::string, Index> const &next_config_id() const;

  std::set<ConfigurationRecord> &data();

  std::set<ConfigurationRecord> const &data() const;

 private:
  std::set<ConfigurationRecord> m_data;

  // map of supercell_name -> next id to assign to a new Configuration
  std::map<std::string, Index> m_next_config_id;
};

/// \brief Make a map for finding ConfigurationRecord by configuration_name
std::map<std::string, ConfigurationRecord const *>
make_index_by_configuration_name(
    std::set<ConfigurationRecord> const &configurations);

}  // namespace config
}  // namespace CASM

#endif
