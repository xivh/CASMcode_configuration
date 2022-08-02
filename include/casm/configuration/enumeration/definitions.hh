// The `casm/configuration/enumeration` module supports
// enumeration of Configuration, IntegralCluster, OccEvents,
// etc.
//
// This module provides:
// - Enumeration methods:
//   - ConfigEnumAllOccupations: for occupation enumeration
// - Filters:
//   - ConfigurationFilter:
//     - AllConfigurationFilter: allow all configurations
//     - UniqueConfigurationFilter: allow primitive, canonical
//       configurations only
//     - GenericConfigurationFilter: customizable filter via
//       std::function and flags to allow only primitive and
//       canonical configurations
//     - ChainedConfigurationFilter: allow only configurations
//       that pass multiple filters
//
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
// - CASMcode_configuration
