#include "casm/configuration/enumeration/perturbations.hh"

#include "casm/configuration/ConfigIsEquivalent.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/enumeration/ConfigEnumAllOccupations.hh"
#include "casm/configuration/enumeration/background_configuration.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/sym_info/definitions.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
namespace config {

/// \brief Make the distinct clusters of sites, taking into account the
///     background configuration symmetry
///
/// \param background, The background
/// \param orbits_as_indices, The orbits in the infinite crystal,
///     converted to linear supercell site indcies
std::set<std::set<Index>> make_distinct_cluster_sites(
    Configuration const &background,
    std::vector<std::set<std::set<Index>>> const &orbits_as_indices) {
  /// Inverse permutations can be used to transform
  /// linear site indices.
  /// Keep only event group operations that also keep
  /// background configuration invariant.
  std::vector<sym_info::Permutation> indices_group_rep;
  ConfigIsEquivalent is_background_invariant(background);
  auto it = SupercellSymOp::begin(background.supercell);
  auto end = SupercellSymOp::end(background.supercell);
  while (it != end) {
    if (is_background_invariant(*it)) {
      indices_group_rep.push_back(sym_info::inverse(it->combined_permute()));
    }
    ++it;
  }
  auto copy_apply_f = [](sym_info::Permutation const &perm,
                         std::set<Index> const &site_indices) {
    std::set<Index> new_site_indices;
    for (Index l : site_indices) {
      new_site_indices.insert(perm[l]);
    }
    return new_site_indices;
  };

  /// Generate new orbit generators.
  /// A generator is the canonical element from an orbit.
  /// These will take into account background configuration and
  /// supercell periodic boundary conditions
  std::set<std::set<Index>> distinct_cluster_sites;
  for (auto const &orbit : orbits_as_indices) {
    std::set<std::set<Index>> tmp = group::make_orbit_generators(
        orbit, indices_group_rep.begin(), indices_group_rep.end(),
        std::less<std::set<Index>>(), copy_apply_f);
    for (auto const &cluster_sites : tmp) {
      distinct_cluster_sites.insert(cluster_sites);
    }
  }

  return distinct_cluster_sites;
}

/// \brief Make configurations that are distinct occupation perturbations
std::set<Configuration> make_distinct_perturbations(
    Configuration const &background,
    std::set<std::set<Index>> const &distinct_cluster_sites) {
  std::set<Configuration> distinct_perturbations;
  auto begin = SupercellSymOp::begin(background.supercell);
  auto end = SupercellSymOp::end(background.supercell);
  for (auto const &cluster_sites : distinct_cluster_sites) {
    ConfigEnumAllOccupations enumerator(background, cluster_sites);
    while (enumerator.is_valid()) {
      distinct_perturbations.emplace(
          make_canonical_form(enumerator.value(), begin, end));
      enumerator.advance();
    }
  }
  return distinct_perturbations;
}

/// \brief Make the distinct clusters of sites, taking into account the
///     event group, supercell, and background configuration symmetry
///
/// \param background, The background
/// \param event_sites Linear sites indices of the cluster of sites that
///     change during the event
/// \param occ_init Initial occupation on sites
/// \param occ_final Final occupation on sites
/// \param event_group The SupercellSymOp consistent with
///     the supercell of the background configuration that leave the
///     event invariant
/// \param local_orbits_as_indices, The local orbits in the infinite crystal,
///     converted to linear supercell site indcies
std::set<std::set<Index>> make_distinct_local_cluster_sites(
    Configuration const &background, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group,
    std::vector<std::set<std::set<Index>>> const &local_orbits_as_indices) {
  /// Inverse permutations can be used to transform
  /// linear site indices.
  /// Keep only event group operations that also keep
  /// background configuration + event combination invariant.

  Configuration config_init = copy_apply_occ(background, event_sites, occ_init);
  Configuration config_final =
      copy_apply_occ(background, event_sites, occ_final);

  std::vector<sym_info::Permutation> indices_group_rep;
  // ConfigIsEquivalent is_background_invariant(background);
  // ConfigIsEquivalent init_is_equiv(config_init);
  // ConfigIsEquivalent final_is_equiv(config_final);
  for (auto const &op : event_group) {
    // Apply completely, then check:
    Configuration A = copy_apply(op, config_init);
    Configuration B = copy_apply(op, config_final);
    if ((A == config_init && B == config_final) ||
        (B == config_init && A == config_final)) {
      indices_group_rep.push_back(sym_info::inverse(op.combined_permute()));
    }

    // if (is_background_invariant(op)) { TODO: check init and final somehow
    //   indices_group_rep.push_back(sym_info::inverse(op.combined_permute()));
    // }

    // // Incremental apply and check:
    // if ( (init_is_equiv(op) && final_is_equiv(op)) || (final_is_equiv(op,
    // config_init) && init_is_equiv(op, config_final)) ) {
    //   indices_group_rep.push_back(sym_info::inverse(op.combined_permute()));
    // }
  }
  auto copy_apply_f = [](sym_info::Permutation const &perm,
                         std::set<Index> const &site_indices) {
    std::set<Index> new_site_indices;
    for (Index l : site_indices) {
      new_site_indices.insert(perm[l]);
    }
    return new_site_indices;
  };

  /// Generate new orbit generators.
  /// A generator is the canonical element from an orbit.
  /// These will take into account background configuration and
  /// supercell periodic boundary conditions
  std::set<std::set<Index>> distinct_local_cluster_sites;
  for (auto const &orbit : local_orbits_as_indices) {
    std::set<std::set<Index>> tmp = group::make_orbit_generators(
        orbit, indices_group_rep.begin(), indices_group_rep.end(),
        std::less<std::set<Index>>(), copy_apply_f);
    for (auto const &local_cluster_sites : tmp) {
      distinct_local_cluster_sites.insert(local_cluster_sites);
    }
  }

  return distinct_local_cluster_sites;
}

/// \brief Make configurations that are distinct local occupation perturbations
///
/// \param background A particular background configuration
/// \param event_sites Linear sites indices of the cluster of sites that
///     change during the event
/// \param occ_init Initial occupation on sites
/// \param occ_final Final occupation on sites
/// \param event_group The SupercellSymOp consistent with
///     the supercell of the background configuration that leave the
///     event invariant
/// \param distinct_local_cluster_sites Symmetrically distinct local-cluster
///     site indices, on which to generate occupation perturbations
std::set<Configuration> make_distinct_local_perturbations(
    Configuration const &background, std::vector<Index> const &event_sites,
    std::vector<int> const &occ_init, std::vector<int> const &occ_final,
    std::vector<SupercellSymOp> const &event_group,
    std::set<std::set<Index>> const &distinct_local_cluster_sites) {
  std::set<Configuration> distinct_local_perturbations;
  for (auto const &local_cluster_sites : distinct_local_cluster_sites) {
    ConfigEnumAllOccupations enumerator(background, local_cluster_sites);
    while (enumerator.is_valid()) {
      distinct_local_perturbations.emplace(make_canonical_form(
          enumerator.value(), event_sites, occ_init, occ_final, event_group));
      enumerator.advance();
    }
  }
  return distinct_local_perturbations;
}

}  // namespace config
}  // namespace CASM
