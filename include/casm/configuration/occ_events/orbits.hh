#ifndef CASM_occ_events_orbits
#define CASM_occ_events_orbits

#include <set>
#include <vector>

#include "casm/configuration/occ_events/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace occ_events {
struct OccEventCounterParameters;

/// \brief Copy OccEvent and apply symmetry operation transformation
OccEvent prim_periodic_occevent_copy_apply(OccEventRep const &rep,
                                           OccEvent occ_event);

/// \brief Make an orbit of OccEvent, with periodic symmetry of a prim
std::set<OccEvent> make_prim_periodic_orbit(
    OccEvent const &orbit_element,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Generate equivalent OccEvent, translated to origin unit cell,
///     in a given order
std::vector<OccEvent> make_prim_periodic_equivalents(
    OccEvent const &prototype, std::vector<Index> const &symop_indices,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make the phenomenal OccEvent, with correct translation
std::vector<OccEvent> make_phenomenal_occevent(
    OccEvent const &prototype,
    std::vector<Index> const &equivalent_generating_op_indices,
    std::vector<xtal::UnitCell> const &phenomenal_generating_translations,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make groups that leave OccEvent orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_occevent_groups(
    std::set<OccEvent> const &orbit,
    std::shared_ptr<SymGroup const> const &symgroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make the group which leaves an OccEvent invariant
std::shared_ptr<SymGroup const> make_occevent_group(
    OccEvent occ_event, std::shared_ptr<SymGroup const> const &symgroup,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make prototypes of distinct orbits of OccEvent, with periodic
///     symmetry of a prim
std::vector<OccEvent> make_prim_periodic_occevent_prototypes(
    std::shared_ptr<OccSystem const> const &system,
    std::vector<clust::IntegralCluster> const &clusters,
    std::vector<OccEventRep> const &occevent_symgroup_rep,
    OccEventCounterParameters const &params,
    std::vector<OccEvent> const &custom_events = {});

/// \brief Make orbits of OccEvent, with periodic symmetry of a prim
std::vector<std::set<OccEvent>> make_prim_periodic_occevent_orbits(
    std::shared_ptr<OccSystem const> const &system,
    std::vector<clust::IntegralCluster> const &clusters,
    std::vector<OccEventRep> const &occevent_symgroup_rep,
    OccEventCounterParameters const &params,
    std::vector<OccEvent> const &custom_events = {});

}  // namespace occ_events
}  // namespace CASM

#endif
