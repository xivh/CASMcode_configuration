#include "casm/configuration/SupercellSet.hh"

#include <map>
#include <set>

#include "casm/configuration/Supercell.hh"
#include "casm/configuration/canonical_form.hh"
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
                              supercell->superlattice.superlattice())),
      is_canonical(config::is_canonical(*supercell)) {
  if (this->is_canonical) {
    this->canonical_supercell_name = this->supercell_name;
  } else {
    auto canonical_supercell = make_canonical_form(*this->supercell);
    this->canonical_supercell_name =
        make_supercell_name(canonical_supercell->superlattice.prim_lattice(),
                            canonical_supercell->superlattice.superlattice());
  }
}

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
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  auto it = find(transformation_matrix_to_super);
  if (it == end()) {
    auto supercell = std::make_shared<Supercell const>(
        m_prim, transformation_matrix_to_super);
    return m_data.emplace(supercell);
  } else {
    return std::make_pair(it, false);
  }
}

/// \brief Insert a canonical supercell by name
///
/// \param supercell_name The name of a canonical supercell
///
/// \returns Returns a pair consisting of an iterator to the inserted element
///     (or to the element that prevented the insertion) and a bool value set to
///     true if and only if the insertion took place.
///
/// Notes:
/// - Throws if `supercell_name` is not the name of the canonical equivalent
/// supercell.
///
std::pair<SupercellSet::iterator, bool> SupercellSet::insert_canonical(
    std::string supercell_name) {
  auto it = find_canonical_by_name(supercell_name);
  if (it == end()) {
    auto supercell = std::make_shared<Supercell const>(
        m_prim, make_superlattice_from_supercell_name(
                    m_prim->basicstructure->lattice(), supercell_name));
    auto canonical_supercell = make_canonical_form(*supercell);
    auto result = m_data.emplace(canonical_supercell);
    if (result.first->canonical_supercell_name != supercell_name) {
      throw std::runtime_error(
          "Error in SupercellSet::insert_canonical: supercell_name is not the "
          "canonical supercell name");
    }
    return result;
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

SupercellSet::const_iterator SupercellSet::find(
    Eigen::Matrix3l const &transformation_matrix_to_super) const {
  auto it = this->begin();
  auto end = this->end();
  for (; it != end; ++it) {
    if (it->supercell->superlattice.transformation_matrix_to_super() ==
        transformation_matrix_to_super) {
      return it;
    }
  }
  return end;
}

SupercellSet::const_iterator SupercellSet::find_canonical_by_name(
    std::string name) const {
  auto it = this->begin();
  auto end = this->end();
  for (; it != end; ++it) {
    if (it->is_canonical && it->supercell_name == name) {
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

SupercellSet::size_type SupercellSet::count(
    Eigen::Matrix3l const &transformation_matrix_to_super) const {
  if (find(transformation_matrix_to_super) != end()) {
    return 1;
  }
  return 0;
}

SupercellSet::size_type SupercellSet::count_canonical_by_name(
    std::string name) const {
  SupercellSet::size_type n = 0;
  auto it = this->begin();
  auto end = this->end();
  for (; it != end; ++it) {
    if (it->is_canonical && it->supercell_name == name) {
      ++n;
    }
  }
  return n;
}

SupercellSet::const_iterator SupercellSet::erase(const_iterator it) {
  return m_data.erase(it);
}

SupercellSet::size_type SupercellSet::erase(
    std::shared_ptr<Supercell const> supercell) {
  return m_data.erase(SupercellRecord(supercell));
}

SupercellSet::size_type SupercellSet::erase(SupercellRecord const &record) {
  return m_data.erase(record);
}

SupercellSet::size_type SupercellSet::erase(
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  auto it = find(transformation_matrix_to_super);
  if (it == end()) {
    return 0;
  }
  m_data.erase(it);
  return 1;
}

SupercellSet::size_type SupercellSet::erase_canonical_by_name(
    std::string name) {
  auto it = find_canonical_by_name(name);
  if (it == end()) {
    return 0;
  }
  m_data.erase(it);
  return 1;
}

std::set<SupercellRecord> &SupercellSet::data() { return m_data; }

std::set<SupercellRecord> const &SupercellSet::data() const { return m_data; }

std::map<std::string, SupercellRecord const *>
make_index_by_canonical_supercell_name(
    std::set<SupercellRecord> const &supercells) {
  std::map<std::string, SupercellRecord const *> result;
  for (auto const &s : supercells) {
    if (s.is_canonical) {
      result.emplace(s.supercell_name, &s);
    }
  }
  return result;
}

/// \brief Find or add canonical supercell by name
SupercellRecord const *find_or_add_canonical_supercell_by_name(
    std::string const &supercell_name, std::set<SupercellRecord> &supercells,
    std::map<std::string, SupercellRecord const *> &index_by_supercell_name,
    std::shared_ptr<Prim const> const &prim) {
  SupercellRecord const *s = nullptr;
  auto it = index_by_supercell_name.find(supercell_name);
  if (it == index_by_supercell_name.end()) {
    auto supercell = std::make_shared<Supercell const>(
        prim, make_superlattice_from_supercell_name(
                  prim->basicstructure->lattice(), supercell_name));
    auto canonical_supercell = make_canonical_form(*supercell);
    s = &*supercells.insert(canonical_supercell).first;
    if (supercell_name != s->supercell_name) {
      throw std::runtime_error(
          "Error in find_or_add_canonical_supercell_by_name: supercell_name "
          "mismatch");
    }
    index_by_supercell_name.emplace(supercell_name, s);
  } else {
    s = it->second;
  }
  return s;
}

}  // namespace config
}  // namespace CASM
