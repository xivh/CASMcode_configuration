#ifndef CASM_config_enum_ConfigEnumAllOccupations
#define CASM_config_enum_ConfigEnumAllOccupations

#include "casm/configuration/Configuration.hh"
#include "casm/container/Counter.hh"

namespace CASM {
namespace config {

/// Enumerate over all possible occupations on particular sites in a
/// Configuration
///
/// Example:
/// \code
/// std::set<Configuration> configurations;
/// Configuration background = ...;
/// std::set<Index> sites = ...;
/// ConfigEnumAllOccupations enumerator(background, sites);
/// std::function<bool (Configuration const &)> filter = ...;
/// while (enumerator.is_valid()) {
///   if (filter(enumerator.value())) {
///     configurations.insert(enumerator.value());
///   }
///   enumerator.advance();
/// }
/// \endcode
///
class ConfigEnumAllOccupations {
 public:
  /// \brief Constructor
  ///
  /// \param background Specifies the background configuration.
  /// \param sites A set of site indices where occupant values are enumerated.
  ///     All other sites in the background configuration maintain the
  ///     original value.
  ///
  ConfigEnumAllOccupations(Configuration const &background,
                           std::set<Index> const &sites);

  /// \brief Get the current Configuration
  Configuration const &value() const;

  /// \brief Generate the next Configuration
  void advance();

  /// \brief Return true if `value` is valid, false if no more valid values
  bool is_valid() const;

 private:
  /// The current configuration
  Configuration m_current;

  /// Site index to enumerate on
  std::set<Index> m_sites;

  /// Counter over allowed occupation indices on sites in m_sites
  Counter<std::vector<int> > m_counter;
};

}  // namespace config
}  // namespace CASM

#endif
