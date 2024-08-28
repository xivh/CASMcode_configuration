#include "casm/configuration/occ_events/orbits.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/group/subgroups.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccEventInvariants.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SymType.hh"

namespace CASM {
namespace occ_events {

/// \brief Copy OccEvent and apply symmetry operation transformation
///
/// \param rep, Symmetry operation representation to be applied
/// \param occ_event, OccEvent to transform
///
/// \return OccEvent, sorted and translated to the origin unit cell
///     after applying the symmetry operation transformation
OccEvent prim_periodic_occevent_copy_apply(OccEventRep const &rep,
                                           OccEvent occ_event) {
  if (!occ_event.size()) {
    return occ_event;
  }
  apply(rep, occ_event);
  clust::IntegralCluster cluster = make_cluster(occ_event);
  occ_event -= cluster[0].unitcell();
  standardize(occ_event);
  return occ_event;
}

/// \brief Make an orbit of OccEvent, with periodic symmetry of a prim
///
/// \param orbit_element One OccEvent in the orbit
/// \param occevent_symgroup_rep Symmetry group representation (as
///     OccEventRep)
///
std::set<OccEvent> make_prim_periodic_orbit(
    OccEvent const &orbit_element,
    std::vector<OccEventRep> const &occevent_symgroup_rep) {
  return group::make_orbit(orbit_element, occevent_symgroup_rep.begin(),
                           occevent_symgroup_rep.end(), std::less<OccEvent>(),
                           prim_periodic_occevent_copy_apply);
}

/// \brief Generate equivalent OccEvent, translated to origin unit cell,
///     in a given order
///
/// \param prototype, A prototype OccEvent
/// \param factor_group_indices, Indices of factor group indices that
///     generate distinct but symmetrically equivalent OccEvent
/// \param lat_column_mat The 3x3 matrix whose columns are the lattice vectors.
/// \param occevent_symgroup_rep Symmetry group representation (as
///     OccEventRep) of the factor group.
///
/// Note:
/// - This method can be used when the order of symmetrically equivalent
///   OccEvent matters, such as for local cluster expansion of OccEvent
///   properties
std::vector<OccEvent> make_prim_periodic_equivalents(
    OccEvent const &prototype, std::vector<Index> const &symop_indices,
    std::vector<OccEventRep> const &occevent_symgroup_rep) {
  std::vector<OccEvent> equivalents;
  for (Index i : symop_indices) {
    equivalents.emplace_back(
        prim_periodic_occevent_copy_apply(occevent_symgroup_rep[i], prototype));
  }
  return equivalents;
}

/// \brief Make the phenomenal OccEvent, with correct translation
///
/// \param prototype The prototype OccEvent
/// \param equivalent_generating_op_indices Factor group operation
///     indices of the symmetry operations used to generate
///     equivalent OccEvent and associated local basis sets
/// \param phenomenal_generating_translations The proper
///     translations (applied after factor group op) for generating
///     OccEvent on the phenomenal cluster sites
/// \param occevent_symgroup_rep Symmetry group representation (as
///     OccEventRep) of the factor group.
///
/// \returns phenomenal_occevent The phenomenal OccEvent, which
///     have the same cluster sites as the phenomenal clusters
///     of the local basis sets, in standardized form.
std::vector<OccEvent> make_phenomenal_occevent(
    OccEvent const &prototype,
    std::vector<Index> const &equivalent_generating_op_indices,
    std::vector<xtal::UnitCell> const &phenomenal_generating_translations,
    std::vector<OccEventRep> const &occevent_symgroup_rep) {
  std::vector<OccEvent> phenomenal_occevent;
  auto const &translations = phenomenal_generating_translations;
  Index i = 0;
  for (Index fg_index : equivalent_generating_op_indices) {
    auto const &occevent_op = occevent_symgroup_rep[fg_index];
    auto const &unitcellcoord_rep = occevent_op.unitcellcoord_rep;
    phenomenal_occevent.push_back(copy_apply(occevent_op, prototype) +
                                  translations[i]);
    ++i;
  }
  for (auto &occ_event : phenomenal_occevent) {
    standardize(occ_event);
  }
  return phenomenal_occevent;
}

/// \brief Make groups that leave OccEvent orbit elements invariant
///
/// \param orbit An OccEvent orbit
/// \param symgroup The symmetry group used to generate the orbit.
/// \param lat_column_mat The 3x3 matrix whose columns are the lattice vectors.
/// \param occevent_symgroup_rep Symmetry group representation (as
///     OccEventRep) of `symgroup`.
///
/// \returns OccEvent invariant groups, where occevent_groups[i] is
///     the SymGroup whose operations leave the sites of the i-th OccEvent in
///     the orbit invariant (up to a permutation/reversal). The head group of
///     the invariant groups is set to be the head group of `symgroup`, which
///     may be `symgroup` itself.
std::vector<std::shared_ptr<SymGroup const>> make_occevent_groups(
    std::set<OccEvent> const &orbit,
    std::shared_ptr<SymGroup const> const &symgroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep) {
  std::shared_ptr<SymGroup const> head_group;
  if (!symgroup->head_group) {
    head_group = symgroup;
  } else {
    head_group = symgroup->head_group;
  }

  // The indices eq_map[i] are the indices of the group
  // elements transform the first element in the orbit into the
  // i-th element in the orbit.
  std::vector<std::vector<Index>> eq_map = group::make_equivalence_map(
      orbit, occevent_symgroup_rep.begin(), occevent_symgroup_rep.end(),
      prim_periodic_occevent_copy_apply);

  // The indices subgroup_indices[i] are the indices of the group
  // elements which leave orbit element i invariant (up to a translation).
  std::vector<group::SubgroupIndices> subgroup_indices =
      group::make_invariant_subgroups(eq_map, *symgroup);

  // The group occevent_groups[i] contains the SymOp corresponding to
  // subgroup_indices[i] and including the translation which keeps
  // the i-th OccEvent invariant
  std::vector<std::shared_ptr<SymGroup const>> occevent_groups;
  auto orbit_it = orbit.begin();
  auto subgroup_indices_it = subgroup_indices.begin();
  auto subgroup_indices_end = subgroup_indices.end();

  while (subgroup_indices_it != subgroup_indices_end) {
    std::vector<xtal::SymOp> occevent_group_elements;
    for (Index j : *subgroup_indices_it) {
      occevent_group_elements.push_back(clust::make_cluster_group_element(
          make_cluster(*orbit_it), lat_column_mat, symgroup->element[j],
          occevent_symgroup_rep[j].unitcellcoord_rep));
    }
    occevent_groups.emplace_back(std::make_shared<SymGroup>(
        head_group, occevent_group_elements, *subgroup_indices_it));
    ++subgroup_indices_it;
    ++orbit_it;
  }
  return occevent_groups;
}

/// \brief Make the group which leaves an OccEvent invariant
///
/// \param occ_event The OccEvent
/// \param symgroup The super group of the occevent group.
/// \param lat_column_mat The 3x3 matrix whose columns are the lattice vectors.
/// \param occevent_symgroup_rep Symmetry group representation (as
///     OccEventRep) of `symgroup`.
///
/// \returns OccEvent invariant group whose operations leave the event
///     invariant (up to a permutation/reversal). The head group of
///     the occevent group is set to be the head group of `symgroup`, which may
///     be `symgroup` itself.
std::shared_ptr<SymGroup const> make_occevent_group(
    OccEvent occ_event, std::shared_ptr<SymGroup const> const &symgroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep) {
  std::shared_ptr<SymGroup const> head_group;
  if (!symgroup->head_group) {
    head_group = symgroup;
  } else {
    head_group = symgroup->head_group;
  }

  if (!occ_event.size()) {
    return std::make_shared<SymGroup>(
        head_group, symgroup->element,
        std::set<Index>(symgroup->head_group_index.begin(),
                        symgroup->head_group_index.end()));
  }

  standardize(occ_event);
  clust::IntegralCluster cluster = make_cluster(occ_event);

  std::vector<xtal::SymOp> elements;
  std::set<Index> indices;
  for (Index i = 0; i < symgroup->element.size(); ++i) {
    OccEvent tocc_event = copy_apply(occevent_symgroup_rep[i], occ_event);
    clust::IntegralCluster tcluster = make_cluster(tocc_event);
    xtal::UnitCell frac_trans = cluster[0].unitcell() - tcluster[0].unitcell();
    tocc_event += frac_trans;
    standardize(tocc_event);

    if (tocc_event == occ_event) {
      xtal::SymOp cart_trans(Eigen::Matrix3d::Identity(),
                             lat_column_mat * frac_trans.cast<double>(), false);
      elements.push_back(cart_trans * symgroup->element[i]);
      indices.insert(i);
    }
  }
  return std::make_shared<SymGroup>(head_group, elements, indices);
}

/// \brief Make prototypes of distinct orbits of OccEvent, with periodic
///     symmetry of a prim
std::vector<OccEvent> make_prim_periodic_occevent_prototypes(
    std::shared_ptr<OccSystem const> const &system,
    std::vector<clust::IntegralCluster> const &clusters,
    std::vector<OccEventRep> const &occevent_symgroup_rep,
    OccEventCounterParameters const &params,
    std::vector<OccEvent> const &custom_events) {
  // function to make an OccEvent canonical
  auto _make_canonical = [&](OccEvent const &event) {
    return group::make_canonical_element(
        event, occevent_symgroup_rep.begin(), occevent_symgroup_rep.end(),
        std::less<OccEvent>(), prim_periodic_occevent_copy_apply);
  };

  typedef std::pair<OccEventInvariants, OccEvent> pair_type;
  CompareOccEvent_f compare_f(system->prim->lattice().tol());
  std::set<pair_type, CompareOccEvent_f> prototype_events(compare_f);

  OccEventCounter counter(system, clusters, params);

  while (!counter.is_finished()) {
    prototype_events.emplace(OccEventInvariants(counter.value(), *system),
                             _make_canonical(counter.value()));
    counter.advance();
  }

  for (auto const &event : custom_events) {
    prototype_events.emplace(OccEventInvariants(event, *system),
                             _make_canonical(event));
  }

  std::vector<OccEvent> result;
  for (auto const &pair : prototype_events) {
    result.emplace_back(pair.second);
  }
  return result;
}

/// \brief Make orbits of OccEvent, with periodic symmetry of a prim
std::vector<std::set<OccEvent>> make_prim_periodic_occevent_orbits(
    std::shared_ptr<OccSystem const> const &system,
    std::vector<clust::IntegralCluster> const &clusters,
    std::vector<OccEventRep> const &occevent_symgroup_rep,
    OccEventCounterParameters const &params,
    std::vector<OccEvent> const &custom_events) {
  std::vector<OccEvent> orbit_prototypes =
      make_prim_periodic_occevent_prototypes(
          system, clusters, occevent_symgroup_rep, params, custom_events);

  // generate orbits from the unique OccEvent
  std::vector<std::set<OccEvent>> orbits;
  for (auto const &event : orbit_prototypes) {
    orbits.emplace_back(make_prim_periodic_orbit(event, occevent_symgroup_rep));
  }

  return orbits;
}

}  // namespace occ_events
}  // namespace CASM
