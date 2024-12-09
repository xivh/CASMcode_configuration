// The `casm/configuration/enumeration` module supports
// enumeration of Configuration
//
// This module provides:
// - Enumeration methods:
//   - ConfigEnumAllOccupations: for occupation enumeration
//   - make_distinct_local_perturbations: for local environment
//     enumeration
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
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_crystallography
// - CASMcode_clexulator
// - CASMcode_configuration (top level)
// - CASMcode_configuration/group
// - CASMcode_configuration/sym_info
// - CASMcode_configuration/clusterography
// - CASMcode_configuration/occ_events
