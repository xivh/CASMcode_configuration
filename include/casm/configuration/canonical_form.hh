#ifndef CASM_config_canonical_form
#define CASM_config_canonical_form

#include <set>

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

struct Configuration;
struct ConfigurationWithProperties;
struct Supercell;
class SupercellSymOp;

// --- Supercell ---

/// \brief Return true if supercell lattice is right-handed lattice in
///    canonical form
bool is_canonical(Supercell const &supercell);

/// \brief Return a shared supercell with right-handed lattice that compares
///     greater to all equivalents with respect to prim point group symmetry
std::shared_ptr<Supercell const> make_canonical_form(
    Supercell const &supercell);

/// \brief Return the supercell with distinct symmetrically equivalent lattices
std::vector<std::shared_ptr<Supercell const>> make_equivalents(
    Supercell const &supercell);

// --- Configuration ---

/// \brief Return true if configuration is in canonical form
template <typename SupercellSymOpIt>
bool is_canonical(Configuration const &configuration, SupercellSymOpIt begin,
                  SupercellSymOpIt end);

/// \brief Return the configuration that compares greater to all equivalents in
///     the same supercell
template <typename SupercellSymOpIt>
Configuration make_canonical_form(Configuration const &configuration,
                                  SupercellSymOpIt begin, SupercellSymOpIt end);

/// \brief Return rep that makes a configuration canonical
template <typename SupercellSymOpIt>
SupercellSymOp to_canonical(Configuration const &configuration,
                            SupercellSymOpIt begin, SupercellSymOpIt end);

/// \brief Return rep that makes a configuration from the canonical
///     configuration
template <typename SupercellSymOpIt>
SupercellSymOp from_canonical(Configuration const &configuration,
                              SupercellSymOpIt begin, SupercellSymOpIt end);

/// \brief Return rep that leave configuration invariant
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    Configuration const &configuration, SupercellSymOpIt begin,
    SupercellSymOpIt end, std::set<std::string> which_dofs = {"all"});

/// \brief Return the distinct symmetrically equivalent configurations (using
///     operations that leave the supercell lattice invariant)
template <typename SupercellSymOpIt>
std::vector<Configuration> make_equivalents(Configuration const &configuration,
                                            SupercellSymOpIt begin,
                                            SupercellSymOpIt end);

/// \brief Return the distinct symmetrically equivalent configurations
///     obtainable by operations consistent with the supercell lattice.
template <typename SupercellSymOpIt>
std::vector<ConfigurationWithProperties> make_equivalents(
    ConfigurationWithProperties const &configuration_with_properties,
    SupercellSymOpIt begin, SupercellSymOpIt end);

/// \brief Find the SupercellSymOp which map equivalent configurations
template <typename SupercellSymOpIt>
std::vector<std::vector<SupercellSymOp>> make_equivalence_map(
    std::vector<Configuration> const &equivalents, SupercellSymOpIt begin,
    SupercellSymOpIt end);

/// \brief Find the SupercellSymOp which map equivalent configurations
template <typename SupercellSymOpIt>
std::vector<std::vector<SupercellSymOp>> make_equivalence_map(
    std::vector<ConfigurationWithProperties> const &equivalents_with_properties,
    SupercellSymOpIt begin, SupercellSymOpIt end);

/// \brief Return true if the operation does not mix given sites and other
/// sites
bool site_indices_are_invariant(SupercellSymOp const &op,
                                std::set<Index> const &site_indices);

/// \brief Return the subgroup of [begin, end] that does not mix given sites and
///     other sites
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    std::set<Index> const &site_indices, SupercellSymOpIt begin,
    SupercellSymOpIt end);

/// \brief Return the subgroup of [begin, end] that leaves configuration
///     invariant and does not mix given sites and other sites
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    Configuration const &configuration, std::set<Index> const &site_indices,
    SupercellSymOpIt begin, SupercellSymOpIt end,
    std::set<std::string> which_dofs = {"all"});

}  // namespace config
}  // namespace CASM

// --- Implementation ---

#include <algorithm>

#include "casm/configuration/ConfigCompare.hh"

namespace CASM {
namespace config {

/// \brief Return true if configuration is in canonical form
///
/// If true, then `configuration` satisfies, for all `rep` in `[begin,
/// end)`:
///     configuration >= copy_apply(rep, configuration)
template <typename SupercellSymOpIt>
bool is_canonical(Configuration const &configuration, SupercellSymOpIt begin,
                  SupercellSymOpIt end) {
  ConfigCompare compare_f(configuration);
  return std::none_of(begin, end, compare_f);
}

/// \brief Return the configuration that compares greater to all equivalents in
///     the same supercell
///
/// The result, `canonical_configuration` satisfies for all `rep` in `[begin,
/// end)`:
///     canonical_configuration >= copy_apply(rep, configuration)
template <typename SupercellSymOpIt>
Configuration make_canonical_form(Configuration const &configuration,
                                  SupercellSymOpIt begin,
                                  SupercellSymOpIt end) {
  return copy_apply(to_canonical(configuration, begin, end), configuration);
}

/// \brief Return rep that makes a configuration canonical
///
/// The result, `rep`, is the first in `[begin, end)` that satisfies:
///     canonical_configuration == copy_apply(rep, configuration)
template <typename SupercellSymOpIt>
SupercellSymOp to_canonical(Configuration const &configuration,
                            SupercellSymOpIt begin, SupercellSymOpIt end) {
  ConfigCompare compare_f(configuration);
  return *std::max_element(begin, end, compare_f);
}

/// \brief Return rep that makes a configuration from the canonical
///     configuration
///
/// The result, `rep`, is the first in `[begin, end)` that satisfies:
///     configuration == copy_apply(rep, canonical_configuration)
template <typename SupercellSymOpIt>
SupercellSymOp from_canonical(Configuration const &configuration,
                              SupercellSymOpIt begin, SupercellSymOpIt end) {
  // simplest version: use the inverse of the first element that results in the
  // canonical form return to_canonical(configuration, begin, end).inverse();

  // alternate version: the lowest index element that transforms canonical form
  // to this
  ConfigCompare compare_f(configuration);
  auto _to_canonical = begin;
  auto _from_canonical = _to_canonical->inverse();
  for (auto it = begin; it < end; ++it) {
    if (compare_f(*_to_canonical, *it)) {
      _to_canonical = it;
      _from_canonical = _to_canonical->inverse();
    }
    // other permutations that result in canonical config may have a lower index
    // inverse
    else if (!compare_f(*it, *_to_canonical)) {
      auto it_inv = it->inverse();
      if (it_inv < _from_canonical) {
        _from_canonical = it_inv;
      }
    }
  }
  return _from_canonical;
}

/// \brief Return rep that leave configuration invariant
///
/// The results, `rep`, are the elements in `[begin, end)` that satisfy:
///     configuration == copy_apply(rep, canonical)
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    Configuration const &configuration, SupercellSymOpIt begin,
    SupercellSymOpIt end, std::set<std::string> which_dofs) {
  std::vector<SupercellSymOp> subgroup;
  ConfigIsEquivalent equal_to_f(configuration, which_dofs);
  std::copy_if(begin, end, std::back_inserter(subgroup), equal_to_f);
  return subgroup;
}

/// \brief Return the distinct symmetrically equivalent configurations
///     obtainable by operations consistent with the supercell lattice.
///
/// The results, `equiv`, are the distinct configuration generated by
/// `equiv = copy_apply(rep, configuration)` for `rep` in `[begin, end)`.
///
/// Note:
/// - This may not generate all configurations that are equivalent as
///   infinite crystals and fit in the supercell, if configuration
///   is not primitive. To generate all equivalents as infinite
///   crystals, use `make_all_super_configurations`.
template <typename SupercellSymOpIt>
std::vector<Configuration> make_equivalents(Configuration const &configuration,
                                            SupercellSymOpIt begin,
                                            SupercellSymOpIt end) {
  std::set<Configuration> equivalents;
  for (auto it = begin; it != end; ++it) {
    equivalents.emplace(copy_apply(*it, configuration));
  }
  return std::vector<Configuration>(equivalents.begin(), equivalents.end());
}

/// \brief Return the distinct symmetrically equivalent configurations
///     obtainable by operations consistent with the supercell lattice.
///
/// The results, `equiv`, are the distinct configuration generated by
/// `equiv = copy_apply(rep, configuration)` for `rep` in `[begin, end)`.
///
/// Note:
/// - Properties are **not** considered in comparisons. It is assumed
///   they are symmetrically consistent with the DoF values.
/// - This may not generate all configurations that are equivalent as
///   infinite crystals and fit in the supercell, if configuration
///   is not primitive. To generate all equivalents as infinite
///   crystals, use `make_all_super_configurations`.
template <typename SupercellSymOpIt>
std::vector<ConfigurationWithProperties> make_equivalents(
    ConfigurationWithProperties const &configuration_with_properties,
    SupercellSymOpIt begin, SupercellSymOpIt end) {
  Configuration const &configuration =
      configuration_with_properties.configuration;
  auto compare = [](std::pair<Configuration, SupercellSymOp> const &A,
                    std::pair<Configuration, SupercellSymOp> const &B) {
    return A.first < B.first;
  };
  std::set<std::pair<Configuration, SupercellSymOp>, decltype(compare)>
      equivalents(compare);
  for (auto it = begin; it != end; ++it) {
    equivalents.emplace(copy_apply(*it, configuration), *it);
  }

  std::vector<ConfigurationWithProperties> equivalents_with_properties;
  if (equivalents.size() == 0) {
    return equivalents_with_properties;
  }
  for (auto const &pair : equivalents) {
    equivalents_with_properties.push_back(
        copy_apply(pair.second, configuration_with_properties));
  }
  return equivalents_with_properties;
}

/// \brief Find the SupercellSymOp which map equivalent configurations
///
/// Generates a table containing which group elements applied to the
/// first equivalent configuration generates each other equivalent.
///
/// \param equivalents The distinct symmetrically equivalent configurations
///     obtainable by operations consistent with the supercell lattice. This
///     may not be all configurations that are equivalent as
///     infinite crystals and fit in the supercell, if the configurations
///     are not primitive.
/// \param begin,end The group used to generate the equivalents
///
/// \returns equivalence_map, The vector equivalence_map[i] is
///     the SupercellSymOp that transform the first element in
///     equivalents into the i-th element in equivalents.
///
template <typename SupercellSymOpIt>
std::vector<std::vector<SupercellSymOp>> make_equivalence_map(
    std::vector<Configuration> const &equivalents, SupercellSymOpIt begin,
    SupercellSymOpIt end) {
  if (equivalents.size() == 0) {
    throw std::runtime_error(
        "Error in make_equivalence_map: equivalents.size() == 0");
  }
  std::vector<std::vector<SupercellSymOp>> equivalence_map;
  equivalence_map.resize(equivalents.size());
  auto equiv_begin = equivalents.begin();
  auto equiv_end = equivalents.end();
  for (auto symop_it = begin; symop_it != end; ++symop_it) {
    auto equiv = copy_apply_f(*symop_it, *equiv_begin);
    Index d = 0;
    auto equiv_it = equiv_begin;
    while (equiv_it != equiv_end && equiv != *equiv_it) {
      ++equiv_it;
      ++d;
    }
    if (equiv_it == equiv_end) {
      throw std::runtime_error("Error in make_equivalence_map: failed");
    }
    equivalence_map[d].push_back(*symop_it);
  }
  return equivalence_map;
}

/// \brief Find the SupercellSymOp which map equivalent configurations
///
/// Generates a table containing which group elements applied to the
/// first equivalent configuration generates each other equivalent.
///
/// \param equivalents_with_properties The distinct symmetrically
///     equivalent configurations obtainable by operations consistent
///     with the supercell lattice. This may not be all configurations
///     that are equivalent as infinite crystals and fit in the
///     supercell, if the configurations are not primitive. Properties
///     are not considered in comparisons; they are assumed to be
///     symmetrically consistent with DoF values.
/// \param begin,end The group used to generate the equivalents
///
/// \returns equivalence_map, The vector equivalence_map[i] is
///     the SupercellSymOp that transform the first element in
///     equivalents into the i-th element in equivalents.
///
template <typename SupercellSymOpIt>
std::vector<std::vector<SupercellSymOp>> make_equivalence_map(
    std::vector<ConfigurationWithProperties> const &equivalents_with_properties,
    SupercellSymOpIt begin, SupercellSymOpIt end) {
  if (equivalents_with_properties.size() == 0) {
    throw std::runtime_error(
        "Error in make_equivalence_map: equivalents_with_properties.size() == "
        "0");
  }
  std::vector<std::vector<SupercellSymOp>> equivalence_map;
  equivalence_map.resize(equivalents_with_properties.size());
  auto equiv_begin = equivalents_with_properties.begin();
  auto equiv_end = equivalents_with_properties.end();
  for (auto symop_it = begin; symop_it != end; ++symop_it) {
    auto equiv = copy_apply_f(*symop_it, *equiv_begin);
    Index d = 0;
    auto equiv_it = equiv_begin;
    while (equiv_it != equiv_end &&
           equiv.configuration != equiv_it->configuration) {
      ++equiv_it;
      ++d;
    }
    if (equiv_it == equiv_end) {
      throw std::runtime_error(
          "Error in make_equivalence_map (with properties): failed");
    }
    equivalence_map[d].push_back(*symop_it);
  }
  return equivalence_map;
}

/// \brief Return the subgroup of [begin, end] that does not mix given sites and
///     other sites
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    std::set<Index> const &site_indices, SupercellSymOpIt begin,
    SupercellSymOpIt end) {
  std::vector<SupercellSymOp> invariant_subgroup;

  for (auto it = begin; it != end; ++it) {
    if (site_indices_are_invariant(*it, site_indices)) {
      invariant_subgroup.push_back(*it);
    }
  }
  return invariant_subgroup;
}

/// \brief Return the subgroup of [begin, end] that leaves configuration
///     invariant and does not mix given sites and other sites
template <typename SupercellSymOpIt>
std::vector<SupercellSymOp> make_invariant_subgroup(
    Configuration const &configuration, std::set<Index> const &site_indices,
    SupercellSymOpIt begin, SupercellSymOpIt end,
    std::set<std::string> which_dofs) {
  std::vector<SupercellSymOp> config_factor_group =
      make_invariant_subgroup(configuration, begin, end, which_dofs);
  return make_invariant_subgroup(site_indices, config_factor_group.begin(),
                                 config_factor_group.end());
}

}  // namespace config
}  // namespace CASM

#endif
