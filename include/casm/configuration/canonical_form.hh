#ifndef CASM_config_canonical_form
#define CASM_config_canonical_form

#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

struct Configuration;
struct Supercell;
class SupercellSymOp;

// --- Supercell ---

/// \brief Return true if supercell lattice is in canonical form
bool is_canonical(Supercell const &supercell);

/// \brief Return a shared supercell that compares greater to all equivalents
///     with respect to prim point group symmetry
std::shared_ptr<Supercell const> make_canonical_form(
    Supercell const &supercell);

/// \brief Return SymOp that makes a supercell lattice canonical
SymOp to_canonical(Supercell const &supercell);

/// \brief Return SymOp that makes a supercell lattice from the canonical
///     supercell lattice
SymOp from_canonical(Supercell const &supercell);

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
    SupercellSymOpIt end);

/// \brief Return the distinct symmetrically equivalent configurations (using
///     operations that leave the supercell lattice invariant)
///
template <typename SupercellSymOpIt>
std::vector<Configuration> make_equivalents(Configuration const &configuration,
                                            SupercellSymOpIt begin,
                                            SupercellSymOpIt end);

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
  auto _from_canonical = _to_canonical.inverse();
  for (auto it = begin; it < end; ++it) {
    if (compare_f(_to_canonical, it)) {
      _to_canonical = it;
      _from_canonical = _to_canonical.inverse();
    }
    // other permutations that result in canonical config may have a lower index
    // inverse
    else if (!compare_f(it, _to_canonical)) {
      auto it_inv = it.inverse();
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
    SupercellSymOpIt end) {
  std::vector<SupercellSymOp> subgroup;
  ConfigIsEquivalent equal_to_f(configuration);
  std::copy_if(begin, end, std::back_inserter(subgroup), equal_to_f);
  return subgroup;
}

/// \brief Return the distinct symmetrically equivalent configurations
///
/// The results, `equiv`, are the distinct configuration generated by
/// `equiv = copy_apply(rep, configuration)` for `rep` in `[begin, end)`.
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

}  // namespace config
}  // namespace CASM

#endif
