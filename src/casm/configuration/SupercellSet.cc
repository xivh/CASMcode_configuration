#include "casm/configuration/SupercellSet.hh"

#include <map>
#include <set>

#include "casm/configuration/Supercell.hh"
#include "casm/configuration/definitions.hh"
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

bool SupercellRecord::operator<(SupercellRecord const &rhs) const {
  return *this->supercell < *rhs.supercell;
}

SupercellSet::SupercellSet(std::shared_ptr<Prim const> const &_prim)
    : m_prim(_prim), m_data() {
  if (m_prim == nullptr) {
    throw std::runtime_error("Error constructing SupercellSet: prim is empty");
  }
}

std::shared_ptr<Prim const> SupercellSet::prim() const { return m_prim; }

bool SupercellSet::empty() const { return m_data.empty(); }

SupercellSet::size_type SupercellSet::size() const { return m_data.size(); }

void SupercellSet::clear() { m_data.clear(); }

SupercellSet::const_iterator SupercellSet::begin() const {
  return m_data.begin();
}

SupercellSet::const_iterator SupercellSet::end() const { return m_data.end(); }

std::pair<SupercellSet::iterator, bool> SupercellSet::insert(
    std::shared_ptr<Supercell const> supercell) {
  return m_data.emplace(supercell);
}

std::pair<SupercellSet::iterator, bool> SupercellSet::insert(
    SupercellRecord const &record) {
  return m_data.insert(record);
}

std::pair<SupercellSet::iterator, bool> SupercellSet::insert(
    std::string supercell_name) {
  auto it = find_by_name(supercell_name);
  if (it == end()) {
    auto supercell = std::make_shared<Supercell const>(
        m_prim, make_superlattice_from_supercell_name(
                    m_prim->basicstructure->lattice(), supercell_name));
    return m_data.emplace(supercell);
  } else {
    return std::make_pair(it, false);
  }
}

SupercellSet::const_iterator SupercellSet::find(
    std::shared_ptr<Supercell const> supercell) const {
  return m_data.find(SupercellRecord(supercell));
}

SupercellSet::const_iterator SupercellSet::find(
    SupercellRecord const &record) const {
  return m_data.find(record);
}

SupercellSet::const_iterator SupercellSet::find_by_name(
    std::string name) const {
  auto it = this->begin();
  auto end = this->end();
  for (; it != end; ++it) {
    if (it->supercell_name == name) {
      return it;
    }
  }
  return end;
}

SupercellSet::size_type SupercellSet::count(
    std::shared_ptr<Supercell const> supercell) const {
  return m_data.count(SupercellRecord(supercell));
}

SupercellSet::size_type SupercellSet::count(
    SupercellRecord const &record) const {
  return m_data.count(record);
}

SupercellSet::size_type SupercellSet::count_by_name(std::string name) const {
  if (find_by_name(name) != end()) {
    return 1;
  }
  return 0;
}

SupercellSet::size_type SupercellSet::erase(
    std::shared_ptr<Supercell const> supercell) {
  return m_data.erase(SupercellRecord(supercell));
}

SupercellSet::size_type SupercellSet::erase(SupercellRecord const &record) {
  return m_data.erase(record);
}

SupercellSet::size_type SupercellSet::erase_by_name(std::string name) {
  auto it = find_by_name(name);
  if (it == end()) {
    return 0;
  }
  m_data.erase(it);
  return 1;
}

std::set<SupercellRecord> &SupercellSet::data() { return m_data; }

std::set<SupercellRecord> const &SupercellSet::data() const { return m_data; }

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
    std::string const &supercell_name, std::set<SupercellRecord> &supercells,
    std::map<std::string, SupercellRecord const *> &index_by_supercell_name,
    std::shared_ptr<Prim const> const &prim) {
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
}  // namespace CASM
