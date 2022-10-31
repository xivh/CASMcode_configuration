// The `casm/configuration/clusterography` module supports
// construction of orbits of clusters
//
// Primarily, this purpose of this module is to provide:
// - IntegralCluster a cluster of UnitCellCoord
// - make_prim_periodic_orbit: function to construct an orbit of
//   IntegralCluster obeying the periodic translation symmetry of
//   a primitive crystal structure
// - make_local_orbit: function to construct an orbit of
//   IntegralCluster without translation symmetry
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_crystallography
// - CASMcode_configuration/group
// - CASMcode_configuration/sym_info

#ifndef CASM_clust_definitions
#define CASM_clust_definitions

#include <functional>
#include <string>

namespace CASM {
template <typename T>
struct traits;

namespace xtal {
class BasicStructure;
class Site;
struct SymOp;
class UnitCell;
class UnitCellCoord;
struct UnitCellCoordRep;
class UnitCellCoordIndexConverter;
}  // namespace xtal

namespace group {
template <typename ElementType>
struct Group;
}

namespace clust {

typedef long Index;
typedef std::string DoFKey;

class ClusterInvariants;
struct ClusterSpecs;
template <typename Base>
class GenericCluster;
class IntegralCluster;
struct IntegralClusterOrbitGenerator;

/// \brief A group::Group of xtal::SymOp
typedef group::Group<xtal::SymOp> SymGroup;

/// A SiteFilterFunction returns true if a Site should be included and false if
/// it should be excluded
typedef std::function<bool(xtal::Site)> SiteFilterFunction;

/// A ClusterFilterFunction returns true if an IntegralCluster should be
/// included and false if it should be excluded
typedef std::function<bool(ClusterInvariants const &, IntegralCluster const &)>
    ClusterFilterFunction;

/// A CandidateSitesFunction generates a vector of UnitCellCoord from a
/// Structure and SiteFilterFuntion
typedef std::function<std::vector<xtal::UnitCellCoord>(
    xtal::BasicStructure const &, SiteFilterFunction)>
    CandidateSitesFunction;

}  // namespace clust
}  // namespace CASM

#endif
