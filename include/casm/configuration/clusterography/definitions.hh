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
