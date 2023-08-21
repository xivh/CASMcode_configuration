# Changelog

All notable changes to `libcasm-configuration` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.0a1] - 2023-08-21

This release creates the libcasm-configuration comparison and enumeration module. It includes:

- Classes for representing supercells, configurations, clusters, and occupation events
- Methods for comparing and enumerating unique configurations, clusters, occupation events, and local environments
- Methods for generating orbits of symmetrically equivalent configurations, clusters, and occupation events
- Methods for copying configurations to make sub- or super-configurations
- Methods for generating symmetry groups, and constructing and applying symmetry representations
- Methods for performing irreducible space decompositions and finding symmetry adapted order parameters

This distribution package libcasm-configuration contains several Python packages of use for configuration comparison and enumeration:

- libcasm.sym_info
- libcasm.irreps
- libcasm.clusterography
- libcasm.configuration
- libcasm.occ_events
- libcasm.enumerate

This package may be installed via pip install, using scikit-build, CMake, and pybind11. This release also includes usage examples and API documentation, built using Sphinx.
