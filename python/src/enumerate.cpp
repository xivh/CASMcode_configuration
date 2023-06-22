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
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/enumeration/MakeOccEventStructures.hh"
#include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/orbits.hh"
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

/// \brief Return (unitcell_index,equivalent_index) of a particular OccEvent
///
/// \param occ_event Input OccEvent, to find the coordinates of
/// \param phenomenal_occevent The phenomenal OccEvent of the equivalent
///     local basis sets
/// \param supercell The supercell in which the OccEvent is located
///
/// \return Returns the coordinate (unitcell_index, equivalent_index) of
///     the input OccEvent. Throws if no match can be found, indicating the
///     input OccEvent is not in the same orbit.
std::tuple<Index, Index> get_occevent_coordinate(
    occ_events::OccEvent occ_event,
    std::vector<occ_events::OccEvent> const &phenomenal_occevent,
    config::Supercell const &supercell) {
  standardize(occ_event);
  clust::IntegralCluster cluster = make_cluster(occ_event);
  cluster.sort();
  for (Index i = 0; i < phenomenal_occevent.size(); ++i) {
    auto phenom_cluster = make_cluster(phenomenal_occevent[i]);
    phenom_cluster.sort();
    xtal::UnitCell trans = phenom_cluster[0].unitcell() - cluster[0].unitcell();
    occ_events::OccEvent translated_occ_event = occ_event + trans;
    if (translated_occ_event == phenomenal_occevent[i]) {
      Index unitcell_index = supercell.unitcell_index_converter(trans);
      Index equivalent_index = i;
      return std::make_tuple(unitcell_index, equivalent_index);
    }
  }
  throw std::runtime_error(
      "Error in get_occevent_coordinate: No match with the phenomenal OccEvent "
      "could be found");
}

std::vector<occ_events::OccEvent> make_phenomenal_occevent(
    occ_events::OccEvent prototype,
    std::vector<clust::IntegralCluster> const &phenomenal_clusters,
    std::vector<Index> const &equivalent_generating_op_indices,
    config::Prim const &prim) {
  clust::IntegralCluster prototype_cluster = make_cluster(prototype);
  prototype_cluster.sort();

  // // get symmetry reps - from xtal::BasicStructure
  // auto factor_group = sym_info::make_factor_group(
  //     xtal_prim);
  // auto unitcellcoord_symgroup_rep =
  //     sym_info::make_unitcellcoord_symgroup_rep(
  //         factor_group->element, xtal_prim);
  // sym_info::OccSymInfo occ_sym_info(
  //     factor_group->element, xtal_prim);
  // auto occevent_symgroup_rep =
  //     occ_events::make_occevent_symgroup_rep(
  //         unitcellcoord_symgroup_rep,
  //         occ_sym_info.occ_symgroup_rep,
  //         occ_sym_info.atom_position_symgroup_rep);

  // get symmetry reps - from config::Prim
  auto const &unitcellcoord_symgroup_rep =
      prim.sym_info.unitcellcoord_symgroup_rep;
  auto occevent_symgroup_rep = occ_events::make_occevent_symgroup_rep(
      unitcellcoord_symgroup_rep, prim.sym_info.occ_symgroup_rep,
      prim.sym_info.atom_position_symgroup_rep);

  // this will throw if there is an inconsistency between prototype,
  // equivalent_generating_op_indices, and phenomenal_clusters,
  // which could happen if the prototype cluster is a different
  // orientation than the one used to generate equivalents
  std::vector<xtal::UnitCell> translations =
      clust::make_phenomenal_generating_translations(
          prototype_cluster, phenomenal_clusters,
          equivalent_generating_op_indices, unitcellcoord_symgroup_rep);

  std::vector<occ_events::OccEvent> phenomenal_occevent =
      occ_events::make_phenomenal_occevent(prototype,
                                           equivalent_generating_op_indices,
                                           translations, occevent_symgroup_rep);

  // sanity check results - this should never throw
  for (Index i = 0; i < phenomenal_clusters.size(); ++i) {
    auto a = phenomenal_clusters[i];
    a.sort();
    auto b = make_cluster(phenomenal_occevent[i]);
    b.sort();
    if (a != b) {
      throw std::runtime_error(
          "Error in libcasm.occ_events.make_phenomenal_occevent: mismatch "
          "between generated phenomenal OccEvent and phenomenal Cluster.");
    }
  }

  return phenomenal_occevent;
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
        "Documented in libcasm.enumerate._methods.py", py::arg("supercell"),
        py::arg("occ_event"), py::arg("motif"), py::arg("local_clusters"));

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

  m.def("make_phenomenal_occevent", &make_phenomenal_occevent,
        R"pbdoc(
      Construct the phenomenal OccEvent for the equivalent local basis sets

      The parameters `phenomenal_clusters` and `equivalent_generating_op_indices`
      can be read from the "equivalents_info.json" file generated when
      the local basis sets are constructed.

      Parameters
      ----------

      prototype : libcasm.occ_events.OccEvent
          The prototype OccEvent used to generate the local basis sets
      phenomenal_clusters : list[libcasm.clusterography.Cluster]
          The phenomenal clusters of the local basis sets, as read
          from "equivalents_info.json"
      equivalent_generating_op_indices : list[int]
          Indices of the factor group operations that generate the
          phenomenal clusters of the local basis sets from the
          prototype, as read from "equivalents_info.json"
      prim : libcasm.config.Prim
          The Prim, used to obtain necessary symmetry info

      Returns
      -------
      phenomenal_occevent : list[libcasm.occ_events.OccEvent]
          The phenomenal OccEvent for the equivalent local basis sets
      )pbdoc",
        py::arg("prototype"), py::arg("phenomenal_clusters"),
        py::arg("equivalent_generating_op_indices"), py::arg("prim"));

  //
  m.def("get_occevent_coordinate", &get_occevent_coordinate,
        R"pbdoc(
      Determine the coordinates `(unitcell_index, equivalent_index)` of a OccEvent

      The coordinates `(unitcell_index, equivalent_index)` are used when
      evaluating local correlations.

      Parameters
      ----------

      occ_event : libcasm.occ_events.OccEvent
          Input OccEvent, to find the coordinates of
      phenomenal_occevent : List[libcasm.occ_events.OccEvent]
          The phenomenal OccEvent for the equivalent local basis sets
      supercell : libcasm.configuration.Supercell
          The supercell in which the OccEvent is located

      Returns
      -------
      (unitcell_index, equivalent_index) : tuple[int, int]
          The coordinates (unitcell_index, equivalent_index) of the
          input OccEvent. Raises if no match can be found, indicating
          the input OccEvent is not in the same orbit.
      )pbdoc",
        py::arg("occ_event"), py::arg("phenomenal_occevent"),
        py::arg("supercell"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
