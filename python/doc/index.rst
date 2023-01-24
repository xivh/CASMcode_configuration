
libcasm-configuration: The Python interface to CASM configuration library
=========================================================================

The libcasm.configuration module is a Python interface to the classes and methods in the CASM::config namespace of the CASM C++ libraries. This includes:

- Classes for representing the primitive crystal structure, symmetry representations, supercells, and configurations
- Methods for applying symmetry to configurations, comparing configurations, and finding the canonical form
- Methods for enumerating unique configurations
- A class for representing clusters of sites
- Methods for generating orbits of symmetrically equivalent clusters
- A class for representing an occupation modifying event for kinetic Monte Carlo calculations
- Methods for finding the symmetry of events and generating orbits of symmetrically equivalent events
- Methods for enumerating unique events
- Methods for enumerating local environments around events


Documentation
=============

.. toctree::
    :maxdepth: 2

    Installation <installation>
    Usage <usage>
    Reference <reference/libcasm/index>
    Bibliography <bibliography>
    About libcasm-configuration <about>


libcasm-configuration is available on GitHub_.

.. _GitHub: https://github.com/prisms-center/CASMcode_configuration
