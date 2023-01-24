#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// nlohmann::json binding
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "pybind11_json/pybind11_json.hpp"

// CASM
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/enumeration/MakeOccEventStructures.hh"
#include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/SimpleStructure.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

std::vector<config::Configuration> make_all_distinct_local_perturbations(
    std::shared_ptr<config::Supercell const> const &supercell,
    occ_events::OccEvent const &occ_event, config::Configuration const &motif,
    std::vector<clust::IntegralCluster> const &local_clusters) {
  std::set<clust::IntegralCluster> _local_clusters(local_clusters.begin(),
                                                   local_clusters.end());

  auto event_prim_info = std::make_shared<config::OccEventPrimInfo const>(
      supercell->prim, occ_event);
  config::OccEventSupercellInfo f(event_prim_info, supercell);
  std::set<config::Configuration> all =
      f.make_all_distinct_local_perturbations(motif, _local_clusters);
  return std::vector<config::Configuration>(all.begin(), all.end());
}

std::vector<xtal::SimpleStructure> make_occevent_simple_structures(
    config::Configuration const &configuration,
    occ_events::OccEvent const &occ_event,
    std::vector<double> interpolation_factors,
    std::shared_ptr<occ_events::OccSystem const> const &system,
    bool skip_event_occupants) {
  config::MakeOccEventStructures make_structure(configuration, occ_event,
                                                system, skip_event_occupants);
  std::vector<xtal::SimpleStructure> results;
  for (double f : interpolation_factors) {
    results.push_back(make_structure(f));
  }
  return results;
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_enumerate, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Methods for enumerating configurations

        libcasm.enumerate
        -----------------

        The libcasm.enumerate package has methods for enumerating configurations.

    )pbdoc";
  py::module::import("libcasm.xtal");
  py::module::import("libcasm.configuration");

  m.def("make_all_distinct_local_perturbations",
        &make_all_distinct_local_perturbations,
        R"pbdoc(
      Construct a subgroup which leaves an event invariant

      Parameters
      ----------
      supercell : libcasm.configuration.Supercell
          The supercell in which local environment configurations will
          be generated.

      occ_event: libcasm.occ_events.OccEvent
          The occupation event.

      motif: libcasm.configuration.Configuration
          The motif configuration is tiled into the supercell to generate
          background configurations. The occ_event is kept in the same position
          while prim factor group symmetry operations are applied to the motif
          configuration to generate all symmetrically distinct combinations of
          the background configuration and occ_event before generating
          perturbations. Only perfect tilings into the supercell are kept.

      local_clusters: list[list[libcasm.clusterography.Cluster]]
          Local clusters, on which the occupation variables will be enumerated
          in order to generate local perturbation configurations. The initial
          cluster in each local-cluster orbit generated using
          :func:`~libcasm.occ_events.make_occevent_cluster_specs` is an
          appropriate input for this parameter. Each local cluster is first used
          to generated local-cluster orbits without regard to the background
          congifuration or supercell. Then, the local-clusters that are distinct
          taking the background configuration and supercell into account are
          perturbed with each possible occupation.

      Returns
      -------
      configurations : libcasm.configuration.Configuration
          The symmetrically distinct perturbations around the
          event on the specified local-clusters, in all distinct
          combinations of the event and motif configuration in
          the chosen supercell.
      )pbdoc",
        py::arg("supercell"), py::arg("occ_event"), py::arg("motif"),
        py::arg("local_clusters"));

  m.def("make_occevent_simple_structures", &make_occevent_simple_structures,
        R"pbdoc(
    Construct Structure along an OccEvent path in a background configuration

    Notes:

    - Only occupation and strain DoF are allowed
    - Vacancies are not included in the resulting strutures
    - Atoms are listed in the same order in all output structures, to
      allow direct input to NEB calculations
    - Raises if multi-atom molecules are allowed in the Prim
    - Raises if event trajectories are not exactly size 2
    - Atom types are sorted according to occ_events::OccSystem::chemical_name_list,
      with atoms involved in the event placed at the beginning of
      each section
    - Atoms are not placed "within" the periodic boundaries if
      the event pathway takes them outside the periodic boundaries

    Parameters
    ----------
    configuration : libcasm.configuration.Configuration
        The background configuration, which sets the occupation on
        all sites not involved in the event.
    occ_event: libcasm.occ_events.OccEvent
        The occupation event.
    interpolation_factors: list[double]
        Interpolation factors, ranging for 0.0 (initial event
        occupation) to 1.0 (final event occupation), specifying
        which structures along the event pathways should be generated.
    system: libcasm.occ_events.OccSystem
        The OccSystem is used to determine output atom type order and
        help with index conversions.
    skip_event_occupants: bool = False
        If True, the occupants involved in the event are not included in
        the output.

    Returns
    -------
    structures : list[xtal.Structure]
        Linearly interpolated structures along the event transformation
        pathway, in the background configuration.
    )pbdoc",
        py::arg("configuration"), py::arg("occ_event"),
        py::arg("interpolation_factors"), py::arg("system"),
        py::arg("skip_event_occupants"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
