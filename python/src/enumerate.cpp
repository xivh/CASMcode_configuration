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
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/configuration/enumeration/ConfigEnumAllOccupations.hh"
#include "casm/configuration/enumeration/MakeOccEventStructures.hh"
#include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/enumeration/background_configuration.hh"
#include "casm/configuration/enumeration/perturbations.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

std::vector<config::Configuration> make_all_distinct_periodic_perturbations(
    std::shared_ptr<config::Supercell const> const &supercell,
    config::Configuration const &motif,
    std::vector<clust::IntegralCluster> const &clusters) {
  auto const &prim = supercell->prim;
  std::vector<std::set<clust::IntegralCluster>> orbits;
  for (auto const &cluster : clusters) {
    orbits.emplace_back(make_prim_periodic_orbit(
        cluster, prim->sym_info.unitcellcoord_symgroup_rep));
  }
  auto orbits_as_indices = clust::make_orbits_as_indices(
      orbits, supercell->unitcellcoord_index_converter);

  std::vector<std::vector<config::Configuration>> super_configurations =
      config::make_all_super_configurations_by_subsets(motif, supercell);

  std::set<config::Configuration> all;
  for (auto const &equiv_configurations : super_configurations) {
    auto const &configuration = equiv_configurations[0];
    auto distinct_cluster_sites =
        config::make_distinct_cluster_sites(configuration, orbits_as_indices);
    std::set<config::Configuration> perturbations =
        config::make_distinct_perturbations(configuration,
                                            distinct_cluster_sites);
    all.insert(perturbations.begin(), perturbations.end());
  }

  return std::vector<config::Configuration>(all.begin(), all.end());
}

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

std::vector<xtal::Lattice> enumerate_superlattices(
    xtal::Lattice const &unit_lattice,
    std::vector<xtal::SymOp> const &point_group, Index max_volume,
    Index min_volume = 1, std::string dirs = std::string("abc"),
    std::optional<Eigen::Matrix3i> unit_cell = std::nullopt,
    bool diagonal_only = false, bool fixed_shape = false) {
  if (!unit_cell.has_value()) {
    unit_cell = Eigen::Matrix3i::Identity();
  }
  xtal::ScelEnumProps enum_props{min_volume,    max_volume + 1,
                                 dirs,          unit_cell.value(),
                                 diagonal_only, fixed_shape};
  xtal::SuperlatticeEnumerator enumerator{unit_lattice, point_group,
                                          enum_props};
  std::vector<xtal::Lattice> superlattices;
  for (auto const &superlat : enumerator) {
    superlattices.push_back(
        xtal::canonical::equivalent(superlat, point_group, unit_lattice.tol()));
  }
  return superlattices;
}

class SuperlatticeEnum {
 public:
  SuperlatticeEnum(xtal::Lattice const &unit_lattice,
                   std::vector<xtal::SymOp> const &point_group,
                   Index max_volume, Index min_volume, std::string dirs,
                   std::optional<Eigen::Matrix3i> unit_cell, bool diagonal_only,
                   bool fixed_shape)
      : m_unit_lattice(unit_lattice),
        m_point_group(point_group),
        m_max_volume(max_volume),
        m_min_volume(min_volume),
        m_dirs(dirs),
        m_unit_cell(unit_cell.has_value() ? unit_cell.value()
                                          : Eigen::Matrix3i::Identity()),
        m_diagonal_only(diagonal_only),
        m_fixed_shape(fixed_shape),
        m_enumerator(
            unit_lattice, point_group,
            xtal::ScelEnumProps{m_min_volume, m_max_volume + 1, m_dirs,
                                m_unit_cell, m_diagonal_only, m_fixed_shape}) {
    m_it = m_enumerator.begin();
    m_end = m_enumerator.end();
    _set_current();
  }

  xtal::Lattice const &value() const { return m_current; }

  void advance() {
    if (m_it != m_end) {
      ++m_it;
    }
    _set_current();
  }

  bool is_valid() const { return m_it != m_end; }

  void reset() {
    m_it = m_enumerator.begin();
    _set_current();
  }

  Index min_volume() const { return m_min_volume; }

  Index max_volume() const { return m_max_volume; }

  std::string dirs() const { return m_dirs; }

  Eigen::Matrix3i unit_cell() const { return m_unit_cell; }

  bool diagonal_only() const { return m_diagonal_only; }

  bool fixed_shape() const { return m_fixed_shape; }

 private:
  void _set_current() {
    if (m_it != m_end) {
      m_current = xtal::canonical::equivalent(*m_it, m_point_group,
                                              m_unit_lattice.tol());
    }
  }

  xtal::Lattice m_unit_lattice;
  std::vector<xtal::SymOp> m_point_group;
  Index m_max_volume;
  Index m_min_volume;
  std::string m_dirs;
  Eigen::Matrix3i m_unit_cell;
  bool m_diagonal_only;
  bool m_fixed_shape;
  xtal::SuperlatticeEnumerator m_enumerator;

  xtal::SuperlatticeIterator m_it;
  xtal::SuperlatticeIterator m_end;
  xtal::Lattice m_current;
};

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
  py::module::import("libcasm.occ_events");

  py::class_<SuperlatticeEnum>(m, "SuperlatticeEnumBase")
      .def(py::init<xtal::Lattice const &, std::vector<xtal::SymOp> const &,
                    Index, Index, std::string, std::optional<Eigen::Matrix3i>,
                    bool, bool>(),
           py::arg("unit_lattice"), py::arg("point_group"),
           py::arg("max_volume"), py::arg("min_volume") = 1,
           py::arg("dirs") = "abc", py::arg("unit_cell") = std::nullopt,
           py::arg("diagonal_only") = false, py::arg("fixed_shape") = false)
      .def("value", &SuperlatticeEnum::value, R"pbdoc(
            Get the current superlattice

            Returns
            -------
            superlattice: xtal.Lattice
                A const reference to the current superlattice.
            )pbdoc",
           py::return_value_policy::reference_internal)
      .def("advance", &SuperlatticeEnum::advance, R"pbdoc(
            Generate the next superlattice
            )pbdoc")
      .def("is_valid", &SuperlatticeEnum::is_valid, R"pbdoc(
          Return True if `value` is valid, False if no more valid values

          Returns
          -------
          is_valid: bool
              True if `value` is valid, False if no more valid values
          )pbdoc")
      .def("reset", &SuperlatticeEnum::reset, R"pbdoc(
            Reset the enumeration to the beginning
            )pbdoc")
      .def("min_volume", &SuperlatticeEnum::min_volume, R"pbdoc(
          Get the minimum volume

          Returns
          -------
          min_volume: int
                  The minimum volume
          )pbdoc")
      .def("max_volume", &SuperlatticeEnum::max_volume, R"pbdoc(
          Get the maximum volume

          Returns
          -------
          max_volume: int
                  The maximum volume
          )pbdoc")
      .def("dirs", &SuperlatticeEnum::dirs, R"pbdoc(
          Get the enumeration directions

          Returns
          -------
          dirs: str
                  The enumeration directions
          )pbdoc")
      .def("unit_cell", &SuperlatticeEnum::unit_cell, R"pbdoc(
          Get the unit cell matrix

          Returns
          -------
          unit_cell: np.ndarray[ np.int[3, 3]]
              The unit cell matrix
          )pbdoc")
      .def("diagonal_only", &SuperlatticeEnum::diagonal_only, R"pbdoc(
          Get the diagonal only flag

          Returns
          -------
          diagonal_only: bool
                  The diagonal only flag
          )pbdoc")
      .def("fixed_shape", &SuperlatticeEnum::fixed_shape, R"pbdoc(
          Get the fixed shape flag

          Returns
          -------
          fixed_shape: bool
                  The fixed shape flag
          )pbdoc");

  py::class_<config::ConfigEnumAllOccupations>(m,
                                               "ConfigEnumAllOccupationsBase")
      .def(py::init<config::Configuration const &, std::set<Index> const &>(),
           py::arg("background"), py::arg("sites"))
      .def("value", &config::ConfigEnumAllOccupations::value, R"pbdoc(
          Get the current Configuration

          Returns
          -------
          config: libcasm.configuration.Configuration
              A const reference to the current Configuration
          )pbdoc",
           py::return_value_policy::reference_internal)
      .def("advance", &config::ConfigEnumAllOccupations::advance, R"pbdoc(
          Generate the next Configuration
          )pbdoc")
      .def("is_valid", &config::ConfigEnumAllOccupations::is_valid, R"pbdoc(
          Return True if `value` is valid, False if no more valid values

          Returns
          -------
          is_valid: bool
              True if `value` is valid, False if no more valid values
          )pbdoc");

  m.def("make_all_distinct_periodic_perturbations",
        &make_all_distinct_periodic_perturbations,
        "Documented in libcasm.enumerate._methods.py", py::arg("supercell"),
        py::arg("motif"), py::arg("clusters"));

  m.def(
      "make_distinct_cluster_sites",
      [](config::Configuration const &background,
         std::vector<std::vector<clust::IntegralCluster>> const &_orbits) {
        // copy vector<vector<IntegralCluster>> -> vector<set<IntegralCluster>>
        std::vector<std::set<clust::IntegralCluster>> orbits;
        for (auto const &_orbit : _orbits) {
          orbits.emplace_back(_orbit.begin(), _orbit.end());
        }

        // convert orbit of IntegralCluster to orbit of linear site indices
        auto supercell = background.supercell;
        xtal::UnitCellCoordIndexConverter const &converter =
            supercell->unitcellcoord_index_converter;
        auto orbits_as_indices =
            clust::make_orbits_as_indices(orbits, converter);

        // find distinct cluster sites in the background configuration
        std::set<std::set<Index>> _distinct_cluster_sites =
            config::make_distinct_cluster_sites(background, orbits_as_indices);

        // copy set<set<Index>> -> vector<set<Index>>
        std::vector<std::set<Index>> distinct_cluster_sites;
        for (auto const &x : _distinct_cluster_sites) {
          distinct_cluster_sites.emplace_back(x);
        }
        return distinct_cluster_sites;
      },
      R"pbdoc(
      Given orbits of clusters in the infinite crystal, make the distinct clusters of
      sites in a particular configuration, as linear site indices

      Parameters
      ----------
      configuration : libcasm.configuration.Configuration
          The background configuration in which distinct clusters will be found.

      orbits: list[list[Cluster]]
          The orbits in the infinite crystal.

      Returns
      -------
      distinct_cluster_sites: list[set[int]]
          The symmetrically distinct clusters in the background configuration,
          where ``distinct_cluster_sites[i]`` is the `i`-th distinct cluster,
          represented as a set of linear site indices in the supercell.
      )pbdoc",
      py::arg("configuration"), py::arg("orbits"));

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

  m.def(
      "make_distinct_local_cluster_sites",
      [](config::Configuration const &configuration,
         occ_events::OccEvent const &occ_event,
         std::vector<config::SupercellSymOp> const &event_group,
         std::vector<std::vector<clust::IntegralCluster>> const
             &_local_orbits) {
        /// \brief Make {cluster, {occ_init, occ_final}} from OccEvent

        auto supercell = configuration.supercell;
        for (auto const &op : event_group) {
          if (op.supercell() != supercell) {
            throw std::runtime_error(
                "Error in libcasm.enumerate.make_distinct_local_cluster_sites: "
                "configuration supercell does not match event_group operation "
                "supercell");
          }
        }

        auto cluster_occupation =
            occ_events::make_cluster_occupation(occ_event);
        std::vector<Index> const &event_sites = to_index_vector(
            cluster_occupation.first, supercell->unitcellcoord_index_converter);
        std::vector<int> const &occ_init = cluster_occupation.second[0];
        std::vector<int> const &occ_final = cluster_occupation.second[1];

        // copy vector<vector<IntegralCluster>> -> vector<set<IntegralCluster>>
        std::vector<std::set<clust::IntegralCluster>> local_orbits;
        for (auto const &_local_orbit : _local_orbits) {
          local_orbits.emplace_back(_local_orbit.begin(), _local_orbit.end());
        }
        // convert to linear site indices
        xtal::UnitCellCoordIndexConverter const &converter =
            supercell->unitcellcoord_index_converter;
        std::vector<std::set<std::set<Index>>> local_orbits_as_indices =
            clust::make_orbits_as_indices(local_orbits, converter);

        std::set<std::set<Index>> _distinct_local_cluster_sites =
            config::make_distinct_local_cluster_sites(
                configuration, event_sites, occ_init, occ_final, event_group,
                local_orbits_as_indices);

        // copy set<set<Index>> -> vector<set<Index>>
        std::vector<std::set<Index>> distinct_local_cluster_sites;
        for (auto const &x : _distinct_local_cluster_sites) {
          distinct_local_cluster_sites.emplace_back(x);
        }
        return distinct_local_cluster_sites;
      },
      R"pbdoc(
      Given orbits of local-clusters in the infinite crystal, make the distinct
      local-clusters around an event in a particular configuration, as linear
      site indices

      Parameters
      ----------
      configuration : libcasm.configuration.Configuration
          The background configuration in which distinct clusters will be found.
      occ_event : libcasm.occ_events.OccEvent
          The event.
      event_group : list[libcasm.configuration.SupercellSymOp]
          The subset of the supercell symmetry group that leaves the event
          invariant.
      local_orbits: list[list[Cluster]]
          The local orbits in the infinite crystal.

      Returns
      -------
      distinct_local_cluster_sites: list[set[int]]
          The symmetrically distinct local-clusters in the background
          configuration, where ``distinct_local_cluster_sites[i]`` is the
          `i`-th distinct cluster, represented as a set of linear site indices
          in the supercell.
      )pbdoc",
      py::arg("configuration"), py::arg("event"), py::arg("event_group"),
      py::arg("local_orbits"));

  m.def(
      "make_distinct_local_perturbations",
      [](config::Configuration const &background,
         occ_events::OccEvent const &occ_event,
         std::vector<config::SupercellSymOp> const &event_group,
         std::vector<std::set<Index>> const &distinct_local_cluster_sites,
         bool allow_subcluster_perturbations) {
        /// \brief Make {cluster, {occ_init, occ_final}} from OccEvent

        auto supercell = background.supercell;
        auto cluster_occupation =
            occ_events::make_cluster_occupation(occ_event);
        std::vector<Index> const &event_sites = to_index_vector(
            cluster_occupation.first, supercell->unitcellcoord_index_converter);
        std::vector<int> const &occ_init = cluster_occupation.second[0];
        std::vector<int> const &occ_final = cluster_occupation.second[1];

        // Choose background + occ_init or background + occ_final as a reference
        config::Configuration reference =
            config::copy_apply_occ(background, event_sites, occ_init);
        config::Configuration config_final =
            config::copy_apply_occ(background, event_sites, occ_final);
        if (config_final > reference) {
          reference = config_final;
        }

        // make distinct local perturbations, with clusters that
        std::set<config::Configuration> distinct_local_perturbations;

        typedef std::tuple<std::vector<Index>, config::Configuration,
                           config::Configuration>
            tuple_type;
        std::vector<tuple_type> results;
        for (auto const &local_cluster_sites : distinct_local_cluster_sites) {
          std::vector<Index> local_cluster_sites_vec(
              local_cluster_sites.begin(), local_cluster_sites.end());
          if (local_cluster_sites.size() == 0) {
            auto result =
                distinct_local_perturbations.emplace(make_canonical_form(
                    background, event_sites, occ_init, occ_final, event_group));
            if (result.second) {
              results.push_back(std::make_tuple(local_cluster_sites_vec,
                                                reference, *result.first));
            }
            continue;
          }
          std::vector<int> initial_occupation;
          std::vector<int> final_occupation;

          // lambda to get occupation on local-cluster sites:
          auto get_occ = [&](config::Configuration const &config,
                             std::vector<int> &occ) {
            occ.resize(local_cluster_sites_vec.size());
            int i = 0;
            for (Index l : local_cluster_sites_vec) {
              occ[i] = config.dof_values.occupation(l);
              ++i;
            }
          };

          get_occ(reference, initial_occupation);

          config::ConfigEnumAllOccupations enumerator(reference,
                                                      local_cluster_sites);
          while (enumerator.is_valid()) {
            get_occ(enumerator.value(), final_occupation);

            bool is_subcluster_perturbation = false;
            for (int i = 0; i < initial_occupation.size(); ++i) {
              if (initial_occupation[i] == final_occupation[i]) {
                is_subcluster_perturbation = false;
                break;
              }
            }
            if (allow_subcluster_perturbations || !is_subcluster_perturbation) {
              auto result = distinct_local_perturbations.emplace(
                  make_canonical_form(enumerator.value(), event_sites, occ_init,
                                      occ_final, event_group));
              if (result.second) {
                results.push_back(std::make_tuple(local_cluster_sites_vec,
                                                  enumerator.value(),
                                                  *result.first));
              }
            }
            enumerator.advance();
          }
        }
        return results;
      },
      R"pbdoc(
      Given orbits of local-clusters in the infinite crystal, make the distinct
      occupation perturbation configurations around an event in a particular
      background configuration.

      Parameters
      ----------
      background : libcasm.configuration.Configuration
          The background configuration in which distinct clusters will be found.
      occ_event : libcasm.occ_events.OccEvent
          The event.
      event_group : list[libcasm.configuration.SupercellSymOp]
          The subset of the supercell symmetry group that leaves the event
          invariant.
      distinct_local_cluster_sites: list[set[int]]
          Local-clusters in the background configuration where perturbations
          should be tried, where ``distinct_local_cluster_sites[i]`` is the
          `i`-th distinct cluster, represented as a set of linear site indices
          in the supercell.

          If `len(distinct_local_cluster_sites[i]) == 0`, the unperturbed
          configuration is included as the only result. Otherwise, only
          perturbations which change the occupation on every site are included
          in the results.
      allow_subcluster_perturbations: bool = False
          If True, output perturbations that only change the occupation on a
          subset of sites. If False (default), only output perturbations that
          change the occupation on all sites.

      Returns
      -------
      results : list[tuple[list[int], libcasm.configuration.Configuration, libcasm.configuration.Configuration]]
          A list of resulting distinct perturbation configurations, as a list
          tuples `(background_index, cluster_sites, config, canonical_config)`,
          where:

          - `cluster_sites`: list[int], The linear site indices for the cluster
            of sites on which the perturbation was applied.
          - `config`: :class:`~libcasm.configuration.Configuration`, The
            configuration with perturbation as applied.
          - `canonical_config`: :class:`~libcasm.configuration.Configuration`,
            The canonical form of the perturbed configuration, as determined by
            application of `event_group`, which is used to find distinct
            configurations.

      )pbdoc",
      py::arg("configuration"), py::arg("event"), py::arg("event_group"),
      py::arg("distinct_local_cluster_sites"),
      py::arg("allow_subcluster_perturbations"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
