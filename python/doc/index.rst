
libcasm-configuration: The Python interface to CASM configuration library
=========================================================================

The libcasm.configuration module is a Python interface to the classes and methods in the CASM C++ libraries that deal with constructing and using configurations, clusters, and occupation events. This includes:


**Configurations:**

- Constructing configurations
- Copying configurations to make sub- or super-configurations
- Enumerating unique configurations


**Clusters:**

- Constructing clusters of sites
- Generating orbits of symmetrically equivalent clusters
- Enumerating unique clusters

**Occupation events:**

- Constructing occupation events for kinetic Monte Carlo calculations
- Enumerating unique events
- Enumerating local environments around events


For each of these objects, there are methods to generate symmetry representations, apply symmetry operations, compare instances and find a canonical form, generate orbits of symmetrically equivalent objects, and enumerate symmetrically unique objects.


Documentation
-------------

.. toctree::
    :maxdepth: 2

    Installation <installation>
    Usage <usage/usage>
    Reference <reference/libcasm/index>
    About <about>


libcasm-configuration is available on GitHub_.

.. _GitHub: https://github.com/prisms-center/CASMcode_configuration
