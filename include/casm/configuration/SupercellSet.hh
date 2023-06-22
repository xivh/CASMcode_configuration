#ifndef CASM_config_SupercellSet
#define CASM_config_SupercellSet

#include <map>
#include <set>

#include "casm/configuration/Supercell.hh"
#include "casm/configuration/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace config {

/// \brief Data structure for holding / reading / writing supercells
struct SupercellRecord : public Comparisons<CRTPBase<SupercellRecord>> {
  SupercellRecord(std::shared_ptr<Supercell const> const &_supercell);

  std::shared_ptr<Supercell const> supercell;

  std::string supercell_name;

  std::string canonical_supercell_name;

  bool is_canonical;

  bool operator<(SupercellRecord const &rhs) const;

 private:
  friend struct Comparisons<CRTPBase<SupercellRecord>>;
};

/// \brief Data structure for holding / reading / writing supercells
class SupercellSet {
 public:
  SupercellSet(std::shared_ptr<Prim const> const &_prim);

  typedef std::set<SupercellRecord>::size_type size_type;
  typedef std::set<SupercellRecord>::iterator iterator;
  typedef std::set<SupercellRecord>::const_iterator const_iterator;

  std::shared_ptr<Prim const> prim() const;

  bool empty() const;

  size_type size() const;

  void clear();

  const_iterator begin() const;

  const_iterator end() const;

  std::pair<iterator, bool> insert(std::shared_ptr<Supercell const> supercell);

  std::pair<iterator, bool> insert(SupercellRecord const &record);

  std::pair<iterator, bool> insert(
      Eigen::Matrix3l const &transformation_matrix_to_super);

  std::pair<iterator, bool> insert_canonical(std::string supercell_name);

  const_iterator find(std::shared_ptr<Supercell const> supercell) const;

  const_iterator find(SupercellRecord const &record) const;

  const_iterator find(
      Eigen::Matrix3l const &transformation_matrix_to_super) const;

  const_iterator find_canonical_by_name(std::string name) const;

  size_type count(std::shared_ptr<Supercell const> supercell) const;

  size_type count(SupercellRecord const &record) const;

  size_type count(Eigen::Matrix3l const &transformation_matrix_to_super) const;

  size_type count_canonical_by_name(std::string name) const;

  const_iterator erase(const_iterator it);

  size_type erase(std::shared_ptr<Supercell const> supercell);

  size_type erase(SupercellRecord const &record);

  size_type erase(Eigen::Matrix3l const &transformation_matrix_to_super);

  size_type erase_canonical_by_name(std::string name);

  std::set<SupercellRecord> &data();

  std::set<SupercellRecord> const &data() const;

 private:
  std::shared_ptr<Prim const> m_prim;
  std::set<SupercellRecord> m_data;
};

/// \brief Make a map for finding canonical SupercellRecord by supercell_name
std::map<std::string, SupercellRecord const *>
make_index_by_canonical_supercell_name(
    std::set<SupercellRecord> const &supercells);

/// \brief Find or add canonical supercell by name
SupercellRecord const *find_or_add_canonical_supercell_by_name(
    std::string const &supercell_name, std::set<SupercellRecord> &supercells,
    std::map<std::string, SupercellRecord const *> &index_by_supercell_name,
    std::shared_ptr<Prim const> const &prim);

}  // namespace config
}  // namespace CASM

#endif
