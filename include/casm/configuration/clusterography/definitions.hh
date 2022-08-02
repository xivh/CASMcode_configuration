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
// - CASMcode_clexulator
// - CASMcode_configuration/group

#ifndef CASM_clust_definitions
#define CASM_clust_definitions

namespace CASM {

namespace xtal {
struct SymOp;
class UnitCell;
class UnitCellCoord;
struct UnitCellCoordRep;
}  // namespace xtal

namespace group {
template <typename ElementType>
struct Group;
}

namespace clust {

typedef long Index;

/// \brief A group::Group of xtal::SymOp
typedef group::Group<xtal::SymOp> SymGroup;

}  // namespace clust
}  // namespace CASM

#endif
