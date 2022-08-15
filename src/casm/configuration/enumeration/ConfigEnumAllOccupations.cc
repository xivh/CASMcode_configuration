#include "casm/configuration/enumeration/ConfigEnumAllOccupations.hh"

namespace CASM {
namespace config {

namespace {  // anonymous

std::vector<int> _make_max_site_occupation(Supercell const &supercell,
                                           std::set<Index> const &sites) {
  auto const &converter = supercell.unitcellcoord_index_converter;
  auto const &basis = supercell.prim->basicstructure->basis();

  auto f = [&](Index site_index) {
    return basis[converter(site_index).sublattice()].occupant_dof().size() - 1;
  };

  std::vector<int> max_site_occupation;
  for (Index i : sites) {
    max_site_occupation.push_back(f(i));
  }

  return max_site_occupation;
}

void _set_occupation(Configuration &configuration, std::set<Index> const &sites,
                     std::vector<int> const &value) {
  Index i = 0;
  for (Index site_index : sites) {
    configuration.dof_values.occupation(site_index) = value[i++];
  }
}

}  // namespace

ConfigEnumAllOccupations::ConfigEnumAllOccupations(
    Configuration const &background, std::set<Index> const &sites)
    : m_current(background),
      m_sites(sites),
      m_counter(std::vector<int>(m_sites.size(), 0),
                _make_max_site_occupation(*m_current.supercell, m_sites),
                std::vector<int>(m_sites.size(), 1)) {
  _set_occupation(m_current, m_sites, m_counter);
}

/// \brief Get the current Configuration
Configuration const &ConfigEnumAllOccupations::value() const {
  return m_current;
}

/// \brief Generate the next Configuration
void ConfigEnumAllOccupations::advance() {
  if (++m_counter) {
    _set_occupation(m_current, m_sites, m_counter);
  }
}

/// \brief Return true if `value` is valid, false if no more values
bool ConfigEnumAllOccupations::is_valid() const { return m_counter.valid(); }

}  // namespace config
}  // namespace CASM
