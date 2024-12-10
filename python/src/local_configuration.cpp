#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
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
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
// #include "casm/configuration/clusterography/orbits.hh"
// #include "casm/configuration/copy_configuration.hh"
// #include "casm/configuration/enumeration/ConfigEnumAllOccupations.hh"
// #include "casm/configuration/enumeration/MakeOccEventStructures.hh"
// #include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/enumeration/background_configuration.hh"
// #include "casm/configuration/enumeration/perturbations.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
// #include "casm/configuration/occ_events/OccSystem.hh"
// #include "casm/configuration/occ_events/orbits.hh"
// #include "casm/crystallography/CanonicalForm.hh"
// #include "casm/crystallography/SimpleStructure.hh"
// #include "casm/crystallography/SuperlatticeEnumerator.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_local_configuration, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Methods for enumerating configurations

        libcasm._local_configuration
        -----------------

        The libcasm._local_configuration package has methods for enumerating
        configurations.

    )pbdoc";
  py::module::import("libcasm.xtal");
  py::module::import("libcasm.configuration");
  py::module::import("libcasm.occ_events");

  m.def(
      "_make_canonical_local_configuration_about_event",
      [](config::Configuration const &configuration,
         occ_events::OccEvent const &occ_event,
         std::vector<config::SupercellSymOp> const &event_group) {
        /// \brief Make {cluster, {occ_init, occ_final}} from OccEvent

        auto supercell = configuration.supercell;
        auto cluster_occupation =
            occ_events::make_cluster_occupation(occ_event);
        std::vector<Index> const &event_sites = to_index_vector(
            cluster_occupation.first, supercell->unitcellcoord_index_converter);
        std::vector<int> const &occ_init = cluster_occupation.second[0];
        std::vector<int> const &occ_final = cluster_occupation.second[1];
        return make_canonical_form(configuration, event_sites, occ_init,
                                   occ_final, event_group);
      },
      R"pbdoc(
      Return the canonical form of a configuration in the context of an OccEvent

      Parameters
      ----------
      local_configuration : libcasm.configuration.Configuration
          The configuration. The occupation on the sites involved in the
          event are not significant - they will be set to either the initial
          or final event occupation.
      occ_event : libcasm.occ_events.OccEvent
          The event.
      event_group : list[libcasm.configuration.SupercellSymOp]
          The subset of the supercell symmetry group that leaves the event
          invariant.

      Returns
      -------
      canonical_local_configuration : libcasm.configuration.Configuration
          The canonical form of the configuration in the context of the event;
          determined by setting the occupation on the sites involved in the
          event to both the initial and final occupation, applying operations
          in the `event_group`, and returning the configuration that compares
          greatest overall.

      )pbdoc",
      py::arg("configuration"), py::arg("event"), py::arg("event_group"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
