#ifndef CASM_config_ConfigIsEquivalent
#define CASM_config_ConfigIsEquivalent

#include "casm/configuration/ConfigDoFIsEquivalent.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/Supercell.hh"

namespace CASM {
namespace config {

/// \brief Class for comparison of Configurations (with the same Supercell)
///
/// - The call operators return the value for equality comparison,
///   and if not equivalent, also store the result for less than comparison
///
class ConfigIsEquivalent {
 public:
  /// Construct with config to be compared against, tolerance for comparison,
  /// and (optional) list of DoFs to compare if _wich_dofs is empty, no dofs
  /// will be compared (default is "all", in which case all DoFs are compared)
  ConfigIsEquivalent(Configuration const &_config, double _tol,
                     std::set<std::string> const &_which_dofs = {"all"});

  ConfigIsEquivalent(Configuration const &_config,
                     std::set<std::string> const &_which_dofs = {"all"});

  Configuration const &config() const;

  /// \brief Returns less than comparison
  ///
  /// - Only valid after call operator returns false
  bool is_less() const;

  /// \brief Check if config == other, store config < other
  ///
  /// - Currently assumes that both Configuration have the same Prim, but may
  ///   have different supercells
  bool operator()(Configuration const &other) const;

  /// \brief Check if config == A*config, store config < A*config
  bool operator()(SupercellSymOp const &A) const;

  /// \brief Check if A*config == B*config, store A*config < B*config
  bool operator()(SupercellSymOp const &A, SupercellSymOp const &B) const;

  /// \brief Check if config == A*other, store config < A*other
  bool operator()(SupercellSymOp const &A, Configuration const &other) const;

  /// \brief Check if A*config == B*other, store A*config < B*other
  bool operator()(SupercellSymOp const &A, SupercellSymOp const &B,
                  Configuration const &other) const;

 private:
  template <typename... Args>
  bool _occupation_is_equivalent(Args &&...args) const;

  Configuration const *m_config;
  Index m_n_sublat;
  bool m_all_dofs;
  bool m_check_occupation;
  bool m_has_aniso_occs;
  Eigen::VectorXi const *m_occupation_ptr;
  std::map<DoFKey, ConfigDoFIsEquivalent::Global> m_global_equivs;
  std::map<DoFKey, ConfigDoFIsEquivalent::Local> m_local_equivs;
  mutable bool m_less;
};

/// Construct with config to be compared against, tolerance for comparison,
/// and (optional) list of DoFs to compare if _wich_dofs is empty, no dofs
/// will be compared (default is "all", in which case all DoFs are compared)
inline ConfigIsEquivalent::ConfigIsEquivalent(
    Configuration const &_config, double _tol,
    std::set<std::string> const &_which_dofs)
    : m_config(&_config),
      m_n_sublat(config().supercell->prim->basicstructure->basis().size()),
      m_all_dofs(_which_dofs.count("all")),
      m_check_occupation(
          (m_all_dofs || _which_dofs.count("occ")) &&
          config().supercell->prim->sym_info.has_occupation_dofs),
      m_has_aniso_occs(config().supercell->prim->sym_info.has_aniso_occs),
      m_occupation_ptr(nullptr) {
  clexulator::ConfigDoFValues const &dof_values = config().dof_values;

  for (auto const &dof : dof_values.global_dof_values) {
    DoFKey const &key = dof.first;
    Eigen::VectorXd const &values = dof.second;
    if (m_all_dofs || _which_dofs.count(key)) {
      m_global_equivs.emplace(std::piecewise_construct,
                              std::forward_as_tuple(key),
                              std::forward_as_tuple(values, key, _tol));
    }
  }

  if (m_check_occupation) {
    m_occupation_ptr = &dof_values.occupation;
  }

  for (auto const &dof : dof_values.local_dof_values) {
    DoFKey const &key = dof.first;
    Eigen::MatrixXd const &values = dof.second;
    if (m_all_dofs || _which_dofs.count(key)) {
      m_local_equivs.emplace(
          std::piecewise_construct, std::forward_as_tuple(key),
          std::forward_as_tuple(values, key, m_n_sublat, _tol));
    }
  }
}

inline ConfigIsEquivalent::ConfigIsEquivalent(
    Configuration const &_config, std::set<std::string> const &_which_dofs)
    : ConfigIsEquivalent(
          _config, _config.supercell->prim->basicstructure->lattice().tol(),
          _which_dofs) {}

inline Configuration const &ConfigIsEquivalent::config() const {
  return *m_config;
}

/// \brief Returns less than comparison
///
/// - Only valid after call operator returns false
inline bool ConfigIsEquivalent::is_less() const { return m_less; }

/// \brief Check if config == other, store config < other
///
/// - Currently assumes that both Configuration have the same Prim, but may
///   have different supercells
inline bool ConfigIsEquivalent::operator()(Configuration const &other) const {
  if (&config() == &other) {
    return true;
  }

  if (config().supercell->prim != other.supercell->prim) {
    throw std::runtime_error(
        "Error comparing Configuration with ConfigIsEquivalent: "
        "Only Configuration with shared prim may be compared this way.");
  }

  clexulator::ConfigDoFValues const &other_dof_values = other.dof_values;

  if (*config().supercell != *other.supercell) {
    m_less = *config().supercell < *other.supercell;
    return false;
  }

  for (auto const &dof_is_equiv_f : m_global_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
    Eigen::VectorXd const &other_values =
        other_dof_values.global_dof_values.at(key);
    if (!f(other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  if (!_occupation_is_equivalent(other_dof_values.occupation)) {
    return false;
  }

  for (auto const &dof_is_equiv_f : m_local_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
    Eigen::MatrixXd const &other_values =
        other_dof_values.local_dof_values.at(key);
    if (!f(other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  return true;
}

/// \brief Check if config == A*config, store config < A*config
inline bool ConfigIsEquivalent::operator()(SupercellSymOp const &A) const {
  for (auto const &dof_is_equiv_f : m_global_equivs) {
    ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
    if (!f(A)) {
      m_less = f.is_less();
      return false;
    }
  }

  if (!_occupation_is_equivalent(A)) {
    return false;
  }

  for (auto const &dof_is_equiv_f : m_local_equivs) {
    ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
    if (!f(A)) {
      m_less = f.is_less();
      return false;
    }
  }

  return true;
}

/// \brief Check if A*config == B*config, store A*config < B*config
inline bool ConfigIsEquivalent::operator()(SupercellSymOp const &A,
                                           SupercellSymOp const &B) const {
  if (A.supercell_factor_group_index() != B.supercell_factor_group_index()) {
    for (auto const &dof_is_equiv_f : m_global_equivs) {
      ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
      if (!f(A, B)) {
        m_less = f.is_less();
        return false;
      }
    }
  }

  if (!_occupation_is_equivalent(A, B)) {
    return false;
  }

  for (auto const &dof_is_equiv_f : m_local_equivs) {
    ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
    if (!f(A, B)) {
      m_less = f.is_less();
      return false;
    }
  }

  return true;
}

/// \brief Check if config == A*other, store config < A*other
inline bool ConfigIsEquivalent::operator()(SupercellSymOp const &A,
                                           Configuration const &other) const {
  clexulator::ConfigDoFValues const &other_dof_values = other.dof_values;

  for (auto const &dof_is_equiv_f : m_global_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
    Eigen::VectorXd const &other_values =
        other_dof_values.global_dof_values.at(key);
    if (!f(A, other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  if (!_occupation_is_equivalent(A, other_dof_values.occupation)) {
    return false;
  }

  for (auto const &dof_is_equiv_f : m_local_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
    Eigen::MatrixXd const &other_values =
        other_dof_values.local_dof_values.at(key);
    if (!f(A, other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  return true;
}

/// \brief Check if A*config == B*other, store A*config < B*other
inline bool ConfigIsEquivalent::operator()(SupercellSymOp const &A,
                                           SupercellSymOp const &B,
                                           Configuration const &other) const {
  clexulator::ConfigDoFValues const &other_dof_values = other.dof_values;

  for (auto const &dof_is_equiv_f : m_global_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Global const &f = dof_is_equiv_f.second;
    Eigen::VectorXd const &other_values =
        other_dof_values.global_dof_values.at(key);
    if (!f(A, B, other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  if (!_occupation_is_equivalent(A, B, other_dof_values.occupation)) {
    return false;
  }

  for (auto const &dof_is_equiv_f : m_local_equivs) {
    DoFKey const &key = dof_is_equiv_f.first;
    ConfigDoFIsEquivalent::Local const &f = dof_is_equiv_f.second;
    Eigen::MatrixXd const &other_values =
        other_dof_values.local_dof_values.at(key);
    if (!f(A, B, other_values)) {
      m_less = f.is_less();
      return false;
    }
  }

  return true;
}

template <typename... Args>
bool ConfigIsEquivalent::_occupation_is_equivalent(Args &&...args) const {
  if (m_check_occupation) {
    if (m_has_aniso_occs) {
      ConfigDoFIsEquivalent::AnisoOccupation f(*m_occupation_ptr, m_n_sublat);
      if (!f(std::forward<Args>(args)...)) {
        m_less = f.is_less();
        return false;
      }
    } else {
      ConfigDoFIsEquivalent::Occupation f(*m_occupation_ptr);
      if (!f(std::forward<Args>(args)...)) {
        m_less = f.is_less();
        return false;
      }
    }
  }
  return true;
}

}  // namespace config
}  // namespace CASM

#endif
