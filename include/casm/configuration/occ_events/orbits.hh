#ifndef CASM_occ_events_orbits
#define CASM_occ_events_orbits

#include <set>
#include <vector>

#include "casm/configuration/occ_events/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace occ_events {

/// \brief Copy OccEvent and apply symmetry operation transformation
OccEvent prim_periodic_occevent_copy_apply(OccEventRep const &rep,
                                           OccEvent occ_event);

/// \brief Find translation that leave OccEvent invariant after
///     transformation, up to a permutation/reversal
xtal::UnitCell prim_periodic_occevent_frac_translation(OccEventRep const &rep,
                                                       OccEvent occ_event);

/// \brief Make an orbit of OccEvent, with periodic symmetry of a prim
std::set<OccEvent> make_prim_periodic_orbit(
    OccEvent const &orbit_element,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make groups that leave OccEvent orbit elements invariant
std::vector<std::shared_ptr<SymGroup const>> make_occevent_groups(
    std::set<OccEvent> const &orbit,
    std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

/// \brief Make the group which leaves an OccEvent invariant
std::shared_ptr<SymGroup const> make_occevent_group(
    OccEvent occ_event, std::shared_ptr<SymGroup const> const &factor_group,
    Eigen::Matrix3d const &lat_column_mat,
    std::vector<OccEventRep> const &occevent_symgroup_rep);

// /// \brief Make orbits of OccEvent, with periodic symmetry of a prim
// std::vector<std::set<OccEvent>> make_prim_periodic_orbits(
//     std::shared_ptr<xtal::BasicStructure const> const &prim,
//     std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
//     sym_info::OccSymGroupRep const &occ_symgroup_rep,
//     sym_info::AtomPositionSymGroupRep const &atom_position_symgroup_rep,
//     SiteFilterFunction site_filter, std::vector<double> const &max_length,
//     std::vector<IntegralClusterOrbitGenerator> const &custom_generators);

}  // namespace occ_events
}  // namespace CASM

#endif
