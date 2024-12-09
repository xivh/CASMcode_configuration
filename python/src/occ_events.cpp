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
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/io/json/OccEventCounter_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
#include "casm/configuration/occ_events/io/stream/OccEventCounter_stream_io.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/LinearIndexConverter.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

std::shared_ptr<occ_events::OccSystem> make_system(
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::optional<std::vector<std::string>> chemical_name_list = std::nullopt,
    std::optional<std::vector<std::string>> vacancy_name_list = std::nullopt) {
  if (!chemical_name_list.has_value()) {
    chemical_name_list =
        occ_events::make_chemical_name_list(*prim, make_factor_group(*prim));
  }

  std::set<std::string> vacancy_name_set = {"Va", "VA", "va"};
  if (vacancy_name_list.has_value()) {
    vacancy_name_set.clear();
    for (auto const &x : vacancy_name_list.value()) {
      vacancy_name_set.insert(x);
    }
  }

  return std::make_shared<occ_events::OccSystem>(prim, *chemical_name_list,
                                                 vacancy_name_set);
}

occ_events::OccEventRep make_occ_event_rep(
    xtal::UnitCellCoordRep const &_unitcellcoord_rep,
    sym_info::OccSymOpRep const &_occupant_rep,
    sym_info::AtomPositionSymOpRep const &_atom_position_rep) {
  return occ_events::OccEventRep(_unitcellcoord_rep, _occupant_rep,
                                 _atom_position_rep);
}

occ_events::OccEvent make_occ_event(
    std::vector<std::vector<occ_events::OccPosition>> trajectories) {
  std::vector<occ_events::OccTrajectory> _trajectories;
  for (auto const &traj : trajectories) {
    _trajectories.emplace_back(traj);
  }
  return occ_events::OccEvent(_trajectories);
}

occ_events::OccPosition make_default_occ_position() {
  return occ_events::OccPosition::molecule(xtal::UnitCellCoord(0, 0, 0, 0), 0);
}

/// \brief Return (unitcell_index,equivalent_index) of a particular OccEvent
///
/// \param occ_event Input OccEvent, to find the coordinates of
/// \param phenomenal_occevent The phenomenal OccEvent of the equivalent
///     local basis sets
/// \param unitcell_index_converter A xtal::UnitCellIndexConverter for the
///     supercell that occ_event is in.
///
/// \return Returns the coordinate (unitcell_index, equivalent_index) of
///     the input OccEvent. Throws if no match can be found, indicating the
///     input OccEvent is not in the same orbit.
std::tuple<Index, Index> get_occevent_coordinate(
    occ_events::OccEvent occ_event,
    std::vector<occ_events::OccEvent> const &phenomenal_occevent,
    xtal::UnitCellIndexConverter const &unitcell_index_converter) {
  standardize(occ_event);
  clust::IntegralCluster cluster = make_cluster(occ_event);
  cluster.sort();
  for (Index i = 0; i < phenomenal_occevent.size(); ++i) {
    auto phenom_cluster = make_cluster(phenomenal_occevent[i]);
    phenom_cluster.sort();
    xtal::UnitCell trans = cluster[0].unitcell() - phenom_cluster[0].unitcell();
    occ_events::OccEvent translated_occ_event = occ_event - trans;
    if (translated_occ_event == phenomenal_occevent[i]) {
      Index unitcell_index = unitcell_index_converter(trans);
      Index equivalent_index = i;
      return std::make_tuple(unitcell_index, equivalent_index);
    }
  }
  throw std::runtime_error(
      "Error in get_occevent_coordinate: No match with the phenomenal OccEvent "
      "could be found");
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_occ_events, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Occupation events, such as diffusive hops or molecular re-orientation

        libcasm.occ_events
        ------------------

        The libcasm.occ_events package contains data structures and methods for specifying and enumerating occupation events, determining their symmetry, and generating orbits.

    )pbdoc";
  py::module::import("libcasm.clusterography");
  py::module::import("libcasm.xtal");
  py::module::import("libcasm.sym_info");

  py::class_<occ_events::OccSystem, std::shared_ptr<occ_events::OccSystem>>(
      m, "OccSystem",
      R"pbdoc(
      Holds index conversion tables used with occupation events
      )pbdoc")
      .def(py::init(&make_system), py::arg("xtal_prim"),
           py::arg("chemical_name_list") = std::nullopt,
           py::arg("vacancy_name_list") = std::nullopt,
           R"pbdoc(

      .. rubric:: Constructor

      Parameters
      ----------
      xtal_prim: libcasm.xtal.Prim
          A :class:`libcasm.xtal.Prim`
      chemical_name_list: Optional[List[str]]=None
          Order of chemical name indices (i.e. :func:`libcasm.xtal.Occupant.name`)
          to use in specifying OccEvents, performing Monte Carlo calculations, etc.
      vacancy_name_list: Optional[List[str]]=None
          Chemical names that should be recognized as vacancies.

      )pbdoc")
      .def(
          "xtal_prim", [](occ_events::OccSystem const &m) { return m.prim; },
          "Return the prim.")
      .def(
          "atom_name_list",
          [](occ_events::OccSystem const &m) { return m.atom_name_list; },
          "Return the atom name list.")
      .def(
          "chemical_name_list",
          [](occ_events::OccSystem const &m) { return m.chemical_name_list; },
          "Return the chemical name list.")
      .def(
          "is_vacancy_list",
          [](occ_events::OccSystem const &m) { return m.is_vacancy_list; },
          "Return a List[bool], where `is_vacancy_list[chemical_name_index]` "
          "indicates if the corresponding chemical is a vacancy.")
      .def(
          "orientation_name_list",
          [](occ_events::OccSystem const &m) {
            return m.orientation_name_list;
          },
          "Names of the unique molecular orientations, as determined from the "
          "keys of :func:`libcasm.xtal.Prim.occupants`.")
      .def_static(
          "from_dict",
          [](const nlohmann::json &data,
             std::shared_ptr<xtal::BasicStructure const> const &prim)
              -> std::shared_ptr<occ_events::OccSystem> {
            // print errors and warnings to sys.stdout
            py::scoped_ostream_redirect redirect;
            jsonParser json{data};
            InputParser<occ_events::OccSystem> event_system_parser(json, prim);
            std::runtime_error error_if_invalid{
                "Error in libcasm.occ_events.OccSystem.from_dict"};
            report_and_throw_if_invalid(event_system_parser, CASM::log(),
                                        error_if_invalid);
            return std::make_shared<occ_events::OccSystem>(
                std::move(*event_system_parser.value));
          },
          R"pbdoc(
          Construct OccSystem from a Python dict

          Parameters
          ----------
          data : dict
              The serialized OccSystem

          xtal_prim : libcasm.xtal.Prim
              A :class:`libcasm.xtal.Prim`

          Returns
          -------
          system : libcasm.occ_events.OccSystem
              The OccSystem
          )pbdoc",
          py::arg("data"), py::arg("xtal_prim"))
      .def(
          "to_dict",
          [](occ_events::OccSystem const &m) -> nlohmann::json {
            jsonParser json;
            to_json(m, json);
            return static_cast<nlohmann::json>(json);
          },
          "Represent the OccSystem as a Python dict.")
      .def("__repr__",
           [](occ_events::OccSystem const &m) {
             std::stringstream ss;
             jsonParser json;
             to_json(m, json);
             ss << json;
             return ss.str();
           })
      .def(
          "get_chemical_name",
          [](occ_events::OccSystem const &self,
             occ_events::OccPosition const &pos) {
            return self.get_chemical_name(pos);
          },
          "Get the chemical name of an occupant", py::arg("pos"))
      .def(
          "get_orientation_name",
          [](occ_events::OccSystem const &self,
             occ_events::OccPosition const &pos) {
            return self.get_orientation_name(pos);
          },

          "Get the orientation name of an occupant", py::arg("pos"))
      .def(
          "get_atom_name",
          [](occ_events::OccSystem const &self,
             occ_events::OccPosition const &pos) {
            return self.get_atom_name(pos);
          },
          "Get the name of an atom component", py::arg("pos"))
      .def("get_chemical_index", &occ_events::OccSystem::get_chemical_index,
           "Get index of occupant in the chemical name list", py::arg("pos"))
      .def("get_orientation_index",
           &occ_events::OccSystem::get_orientation_index,
           "Get index of occupant in the orientation name list", py::arg("pos"))
      .def("get_cartesian_coordinate",
           &occ_events::OccSystem::get_cartesian_coordinate,
           "Get the Cartesian coordinate of an occupant position",
           py::arg("pos"))
      .def("is_indivisible", &occ_events::OccSystem::is_indivisible,
           "Return True if occupant is indivisible, False otherwise",
           py::arg("pos"))
      .def("is_vacancy", &occ_events::OccSystem::is_vacancy,
           "Return True if occupant is a vacancy, False otherwise",
           py::arg("pos"))
      .def("get_occupant", &occ_events::OccSystem::get_occupant,
           "Return the :class:`~libcasm.xtal.Occupant` indicated",
           py::arg("pos"));

  py::class_<occ_events::OccPosition>(m, "OccPosition",
                                      R"pbdoc(
      An atom or molecule position
      )pbdoc")
      .def(py::init(&make_default_occ_position),
           R"pbdoc(

      .. rubric:: Constructor

      Construct a default OccPosition, equal to indicating the first occupant in the
      first basis site in the origin unit cell.

      )pbdoc")
      .def_static(
          "molecule",
          [](xtal::UnitCellCoord const &integral_site_coordinate,
             Index occupant_index) {
            return occ_events::OccPosition::molecule(integral_site_coordinate,
                                                     occupant_index);
          },
          R"pbdoc(
          Construct an OccPosition representing an entire molecule, whether
          single or multi-atom.

          This is equivalent to :func:`~libcasm.occ_events.occupant`.
          )pbdoc",
          py::arg("integral_site_coordinate"), py::arg("occupant_index"))
      .def_static(
          "occupant",
          [](xtal::UnitCellCoord const &integral_site_coordinate,
             Index occupant_index) {
            return occ_events::OccPosition::molecule(integral_site_coordinate,
                                                     occupant_index);
          },
          R"pbdoc(
          Construct an OccPosition representing the entire occupant.

          This is equivalent to :func:`~libcasm.occ_events.molecule`.
          )pbdoc",
          py::arg("integral_site_coordinate"), py::arg("occupant_index"))
      .def_static(
          "atom_component",
          [](xtal::UnitCellCoord const &integral_site_coordinate,
             Index occupant_index, Index atom_position_index) {
            return occ_events::OccPosition::atom(
                integral_site_coordinate, occupant_index, atom_position_index);
          },
          "Construct an OccPosition representing an atom component of a "
          "multi-atom molecule.",
          py::arg("integral_site_coordinate"), py::arg("occupant_index"),
          py::arg("atom_position_index"))
      .def(
          "is_in_reservoir",
          [](occ_events::OccPosition const &pos) {
            return pos.is_in_reservoir;
          },
          "If true, indicates molecule/atom in reservoir. If false, indicates "
          "a "
          "molecule/atom on integral_site_coordinate.")
      .def(
          "is_atom",
          [](occ_events::OccPosition const &pos) { return pos.is_atom; },
          "If true, indicates this tracks an atom component. If false, then "
          "this tracks a molecule position.")
      .def(
          "integral_site_coordinate",
          [](occ_events::OccPosition const &pos) {
            return pos.integral_site_coordinate;
          },
          "If is_in_reservoir() is False: Integral coordinates of site "
          "containing occupant; otherwise invalid.")
      .def(
          "occupant_index",
          [](occ_events::OccPosition const &pos) { return pos.occupant_index; },
          "If is_in_reservoir() is False: Index of occupant in "
          ":func:`~libcasm.xtal.Prim.occ_dof` for sublattice specified by "
          "`integral_site_coordinate`.  If is_in_reservoir() is True: Index "
          "into :func:`~libcasm.occ_events.OccSystem.chemical_name_list` of a "
          "molecule in the reservoir.")
      .def(
          "atom_position_index",
          [](occ_events::OccPosition const &pos) { return pos.occupant_index; },
          "If is_atom() is True and is_in_reservoir() is False: Index of atom "
          "position in the indicated occupant molecule.")
      .def(
          "copy",
          [](occ_events::OccPosition const &self) {
            return occ_events::OccPosition(self);
          },
          R"pbdoc(
          Returns a copy of the OccPosition.
          )pbdoc")
      .def("__copy__",
           [](occ_events::OccPosition const &self) {
             return occ_events::OccPosition(self);
           })
      .def("__deepcopy__",
           [](occ_events::OccPosition const &self, py::dict) {
             return occ_events::OccPosition(self);
           })
      .def_static(
          "from_dict",
          [](const nlohmann::json &data,
             occ_events::OccSystem const &system) -> occ_events::OccPosition {
            jsonParser json{data};
            return jsonConstructor<occ_events::OccPosition>::from_json(json,
                                                                       system);
          },
          R"pbdoc(
          Construct an OccPosition from a Python dict

          Parameters
          ----------
          data : dict
              The serialized OccPosition

          system : libcasm.occ_events.OccSystem
              A :class:`~libcasm.occ_events.OccSystem`

          Returns
          -------
          event : libcasm.occ_events.OccPosition
              The OccPosition
          )pbdoc",
          py::arg("data"), py::arg("system"))
      .def(
          "to_dict",
          [](occ_events::OccPosition const &pos,
             occ_events::OccSystem const &system) -> nlohmann::json {
            jsonParser json;
            to_json(pos, json, system);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the OccPosition as a Python dict

          Parameters
          ----------
          system : libcasm.occ_events.OccSystem
              A :class:`~libcasm.occ_events.OccSystem`

          Returns
          -------
          data : dict
              The OccPosition as a Python dict
          )pbdoc",
          py::arg("system"))
      .def("__repr__", [](occ_events::OccPosition const &pos) {
        std::stringstream ss;
        jsonParser json;
        to_json(pos, json, std::nullopt);
        ss << json;
        return ss.str();
      });

  py::class_<occ_events::OccEventRep> pyOccEventRep(m, "OccEventRep",
                                                    R"pbdoc(
      Symmetry representation for transformating an OccEvent.
    )pbdoc");

  py::class_<occ_events::OccEvent> pyOccEvent(m, "OccEvent",
                                              R"pbdoc(
      OccEvent represents an occupation event, for example the change in
      occupation due to a diffusive hop or molecular re-orientation. The
      occupation change is represented by occupant trajectories.

      Example, 1NN A-Va exchange in an FCC prim:

      .. code-block:: Python

          import libcasm.xtal as xtal
          from libcasm.xtal.prims import FCC as FCC_prim
          from libcasm.occ_events import OccPosition, OccEvent

          prim = FCC_prim(r=1.0, occ_dof=["A", "B", "Va"])

          site1 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[0, 0, 0])
          site2 = xtal.IntegralSiteCoordinate(sublattice=0, unitcell=[1, 0, 0])

          A_occ_index = 0
          Va_occ_index = 2

          A_initial_pos = OccPosition.molecule(site1, A_occ_index)
          A_final_pos = OccPosition.molecule(site2, A_occ_index)
          Va_initial_pos = OccPosition.molecule(site2, Va_occ_index)
          Va_final_pos = OccPosition.molecule(site1, Va_occ_index)

          occ_event = OccEvent([
              [A_initial_pos, A_final_pos],
              [Va_initial_pos, Va_final_pos]
          ])

    )pbdoc");

  pyOccEventRep
      .def(py::init(&make_occ_event_rep),
           R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          integral_site_coordinate_rep: libcasm.xtal.IntegralSiteCoordinateRep
              Symmetry representation for transforming IntegralSiteCoordinate

          occupant_rep: List[List[int]]
              Permutations describe occupant index transformation under symmetry.
              Usage:

                  occupant_index_after =
                      occupant_rep[sublattice_index_before][occupant_index_before]

          atom_position_rep: List[List[List[int]]]
              Permutations describe atom position index transformation under
              symmetry.

              Usage:

                  atom_position_index_after =
                      atom_position_rep[sublattice_index_before][occupant_index_before][atom_position_index_before]

           )pbdoc",
           py::arg("integral_site_coordinate_rep"), py::arg("occupant_rep"),
           py::arg("atom_position_rep"))
      .def(
          "__mul__",
          [](occ_events::OccEventRep const &self,
             occ_events::OccEvent const &event) {
            return copy_apply(self, event);
          },
          py::arg("event"),
          "Creates a copy of the OccEvent and applies the symmetry operation "
          "represented by this OccEventRep.")
      .def(
          "copy",
          [](occ_events::OccEventRep const &self) {
            return occ_events::OccEventRep(self);
          },
          R"pbdoc(
          Returns a copy of the OccEventRep.
          )pbdoc")
      .def("__copy__",
           [](occ_events::OccEventRep const &self) {
             return occ_events::OccEventRep(self);
           })
      .def("__deepcopy__",
           [](occ_events::OccEventRep const &self, py::dict) {
             return occ_events::OccEventRep(self);
           })
      .def("__repr__", [](occ_events::OccEventRep const &self) {
        std::stringstream ss;
        jsonParser json;
        json["site_rep"]["sublattice_after"] =
            self.unitcellcoord_rep.sublattice_index;
        json["site_rep"]["matrix_frac"] = self.unitcellcoord_rep.point_matrix;
        json["site_rep"]["tau_frac"] = jsonParser::array();
        for (auto const &tau_frac : self.unitcellcoord_rep.unitcell_indices) {
          jsonParser tjson;
          to_json(tau_frac, tjson, jsonParser::as_array());
          json["site_rep"]["tau_frac"].push_back(tjson);
        }
        json["occupant_rep"] = self.occupant_rep;
        json["atom_position_rep"] = self.atom_position_rep;
        ss << json;
        return ss.str();
      });

  m.def(
      "make_occevent_symgroup_rep",
      [](std::vector<xtal::SymOp> const &group_elements,
         xtal::BasicStructure const &xtal_prim) {
        return occ_events::make_occevent_symgroup_rep(group_elements,
                                                      xtal_prim);
      },
      R"pbdoc(
      Construct a group representation of OccEventRep

      Parameters
      ----------
      group_elements: List[libcasm.xtal.SymOp]
          Symmetry group elements

      xtal_prim: libcasm.xtal.Prim
          The Prim structure

      Returns
      -------
      occevent_symgroup_rep: List[OccEventRep]
          Group representation for transforming OccEvent
      )pbdoc",
      py::arg("group_elements"), py::arg("xtal_prim"));

  m.def(
      "make_occevent_symgroup_rep_from_existing",
      [](std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep,
         std::vector<sym_info::OccSymOpRep> const &occ_symgroup_rep,
         std::vector<sym_info::AtomPositionSymOpRep> const
             &atom_position_symgroup_rep) {
        return occ_events::make_occevent_symgroup_rep(
            unitcellcoord_symgroup_rep, occ_symgroup_rep,
            atom_position_symgroup_rep);
      },
      R"pbdoc(
      Construct a group representation of OccEventRep

      Parameters
      ----------
      unitcellcoord_symgroup_rep: List[libcasm.xtal.IntegralSiteCoordinateRep]
          The symmetry group representation that describes how IntegralSiteCoordinate transform under symmetry.
      occ_symgroup_rep: List[List[List[int]]]
          Permutations describe occupant index transformation under symmetry. Indices are: group_element_index, sublattice_index_before, and occupant_index_before; and the resulting value is occupant_index_after.
      atom_position_symgroup_rep: List[List[List[List[int]]]]
          Permutations describe atom position index transformation under symmetry. Indices are: group_element_index, sublattice_index_before, occupant_index_before, atom_position_index_before; and the resulting value is atom_position_index_after.

      Returns
      -------
      occevent_symgroup_rep: List[OccEventRep]
          Group representation for transforming OccEvent
      )pbdoc",
      py::arg("unitcellcoord_symgroup_rep"), py::arg("occ_symgroup_rep"),
      py::arg("atom_position_symgroup_rep"));

  pyOccEvent
      .def(py::init(&make_occ_event),
           R"pbdoc(

      .. rubric:: Special Methods

      The multiplication operator ``X = lhs * rhs`` can be used to apply
      :class:`libcasm.xtal.OccEventRep` to OccEvent:

      - ``X=OccEvent``, ``lhs=OccEventRep``, ``rhs=OccEvent``: Copy and
        transform the trajectories, returning the OccEvent of transformed
        trajectories. This transforms the sites involved, and also the
        occupation indices (for anisotropic occupants, or for sublattices with
        different allowed occupation lists) and atom positions (for molecular
        occupants).

      Translate an OccEvent using operators ``+``, ``-``, ``+=``, ``-=``:

      .. code-block:: Python

          import numpy as np
          from libcasm.occ_events import OccEvent

          # event: OccEvent
          translation = np.array([0, 0, 1])

          # translate via `+=`:
          event += translation

          # translate via `-=`:
          event -= translation

          # copy & translate via `+`:
          translated_event = event + translation

          # copy & translate via `-`:
          translated_event = event - translation


      Additional methods:

      - Sort OccEvent by lexicographical order of trajectories using ``<``, ``<=``,
        ``>``, ``>=``, and compare using ``==`` and ``!=``
      - ``for site in cluster``: Iterate over sites
        (:class:`~libcasm.xtal.IntegralSiteCoordinate`) in the cluster.
      - ``if site in cluster``: Check if a cluster contains a site.
      - ``len(cluster)``: Get the number of sites in a cluster.
      - ``site = cluster[i]``: Get the `i`-th site in a cluster (indices
        start at 0).
      - OccEvent may be copied with
        :func:`OccEvent.copy <libcasm.occ_events.OccEvent.copy>`,
        `copy.copy`, or `copy.deepcopy`.

      .. rubric:: Constructor

      Parameters
      ----------
      trajectories: List[List[OccPosition]]=[]
          The occupant trajectories. Usage: `trajectories[i_occupant][0]` is the initial position of the i-th occupant, and `trajectories[i_occupant][1]` is the final position of the i-th occupant. Most methods currently support trajectories of length 2 only (an initial position and a final position).

      )pbdoc",
           py::arg("trajectories") =
               std::vector<std::vector<occ_events::OccPosition>>{})
      .def(
          "size",
          [](occ_events::OccEvent const &event) { return event.size(); },
          "The number of trajectories.")
      .def(
          "trajectories",
          [](occ_events::OccEvent const &event) {
            std::vector<std::vector<occ_events::OccPosition>> trajectories(
                event.size());
            for (Index i_traj = 0; i_traj < event.size(); ++i_traj) {
              for (auto const &pos : event[i_traj].position) {
                trajectories[i_traj].push_back(pos);
              }
            }
            return trajectories;
          },
          R"pbdoc(
          Return the event trajectories

          Returns
          -------
          trajectories: List[List[OccPosition]]=[]
             The occupant trajectories. Usage: `trajectories[i_occupant][0]` is the initial position of the i-th occupant, and `trajectories[i_occupant][1]` is the final position of the i-th occupant. Most methods currently support trajectories of length 2 only (an initial position and a final position).

          )pbdoc")
      .def(
          "__add__",
          [](occ_events::OccEvent const &event,
             Eigen::Vector3l const &translation) {
            return event + translation;
          },
          "Translate the OccEvent by adding unit cell indices")
      .def(
          "__iadd__",
          [](occ_events::OccEvent &event, Eigen::Vector3l const &translation) {
            return event += translation;
          },
          "Translate the OccEvent by adding unit cell indices")
      .def(
          "__sub__",
          [](occ_events::OccEvent const &event,
             Eigen::Vector3l const &translation) {
            return event - translation;
          },
          "Translate the OccEvent by subtracting unit cell indices")
      .def(
          "__isub__",
          [](occ_events::OccEvent &event, Eigen::Vector3l const &translation) {
            return event -= translation;
          },
          "Translate the OccEvent by subtracting unit cell indices")
      .def(
          "sort", [](occ_events::OccEvent &event) { return sort(event); },
          "Sort event trajectories")
      .def(
          "copy_sort",
          [](occ_events::OccEvent const &event) { return copy_sort(event); },
          "Return a copy of the event with sorted trajectories")
      .def(
          "reverse", [](occ_events::OccEvent &event) { return reverse(event); },
          "Reverse event trajectories")
      .def(
          "copy_reverse",
          [](occ_events::OccEvent const &event) { return copy_reverse(event); },
          "Return a copy of the event with reversed trajectories")
      .def(
          "standardize",
          [](occ_events::OccEvent &event) { return standardize(event); },
          R"pbdoc(
          Put event into standardized form with regard to permutation/reversal.

          This is equivalent to:

          .. code-block:: Python

              forward_occ_event = self.copy().sort();
              reverse_occ_event = self.copy().reverse().sort();
              if (reverse_occ_event < occ_event) {
                self = reverse_occ_event;
              }
              else {
                self = forward_occ_event;
              }

          )pbdoc"
          "")
      .def(
          "cluster",
          [](occ_events::OccEvent const &event) {
            return make_cluster_occupation(event).first;
          },
          "The cluster of sites involved in the OccEvent.")
      .def(
          "initial_occupation",
          [](occ_events::OccEvent const &event) {
            return make_cluster_occupation(event).second[0];
          },
          "Occupant indices on each site in the cluster, in the initial "
          "positions. Order of sites is consistent with self.cluster().")
      .def(
          "final_occupation",
          [](occ_events::OccEvent const &event) {
            return make_cluster_occupation(event).second[1];
          },
          "Occupant indices on each site in the cluster, in the final "
          "positions. Order of sites is consistent with self.cluster().")
      .def(py::self < py::self,
           "Compares via lexicographical order of trajectories")
      .def(py::self <= py::self,
           "Compares via lexicographical order of trajectories")
      .def(py::self > py::self,
           "Compares via lexicographical order of trajectories")
      .def(py::self >= py::self, "Compares via lexicographical order of sites")
      .def(py::self == py::self, "True if events are equal")
      .def(py::self != py::self, "True if events are not equal")
      .def(
          "copy",
          [](occ_events::OccEvent const &self) {
            return occ_events::OccEvent(self);
          },
          R"pbdoc(
          Returns a copy of the OccEvent.
          )pbdoc")
      .def("__copy__",
           [](occ_events::OccEvent const &self) {
             return occ_events::OccEvent(self);
           })
      .def("__deepcopy__", [](occ_events::OccEvent const &self,
                              py::dict) { return occ_events::OccEvent(self); })
      .def_static(
          "from_dict",
          [](const nlohmann::json &data,
             occ_events::OccSystem const &system) -> occ_events::OccEvent {
            jsonParser json{data};
            return jsonConstructor<occ_events::OccEvent>::from_json(json,
                                                                    system);
          },
          R"pbdoc(
          Construct an OccEvent from a Python dict

          Parameters
          ----------
          data : dict
              The serialized OccEvent

          system : libcasm.occ_events.OccSystem
              A :class:`~libcasm.occ_events.OccSystem`

          Returns
          -------
          event : libcasm.occ_events.OccEvent
              The OccEvent
          )pbdoc",
          py::arg("data"), py::arg("system"))
      .def(
          "to_dict",
          [](occ_events::OccEvent const &event,
             occ_events::OccSystem const &system, bool include_cluster,
             bool include_cluster_occupation,
             bool include_event_invariants) -> nlohmann::json {
            occ_events::OccEventOutputOptions opt;
            opt.include_cluster = include_cluster;
            opt.include_cluster_occupation = include_cluster_occupation;
            opt.include_event_invariants = include_event_invariants;
            jsonParser json;
            to_json(event, json, system, opt);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the OccEvent as a Python dict

          Parameters
          ----------
          system : libcasm.occ_events.OccSystem
              A :class:`~libcasm.occ_events.OccSystem`

          include_cluster: bool = True
              If True, also include the cluster sites

          include_cluster_occupation: bool = True
              If True, also include the initial and final cluster occupation

          include_event_invariants: bool = True
              If True, also include event invariants: number of trajectories,
              number of each occupant type, and site distances

          Returns
          -------
          data : dict
              The OccEvent as a Python dict
          )pbdoc",
          py::arg("system"), py::arg("include_cluster") = true,
          py::arg("include_cluster_occupation") = true,
          py::arg("include_event_invariants") = true)
      .def("__repr__", [](occ_events::OccEvent const &event) {
        std::stringstream ss;
        occ_events::OccEventOutputOptions opt;
        opt.include_cluster = false;
        opt.include_cluster_occupation = true;
        opt.include_event_invariants = false;
        jsonParser json;
        to_json(event, json, std::nullopt, opt);
        ss << json;
        return ss.str();
      });

  m.def(
      "make_prim_periodic_orbit",
      [](occ_events::OccEvent const &orbit_element,
         std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep) {
        std::set<occ_events::OccEvent> orbit =
            make_prim_periodic_orbit(orbit_element, occevent_symgroup_rep);
        return std::vector<occ_events::OccEvent>(orbit.begin(), orbit.end());
      },
      R"pbdoc(
      Construct an orbit of OccEvent

      The orbit of OccEvent is all distinct OccEvent that are equivalent
      under the provided symmetry group, including one element for all
      OccEvent that are equivalent according to prim translational symmetry.

      Parameters
      ----------
      orbit_element : OccEvent
          One OccEvent in the orbit

      occevent_symgroup_rep: List[OccEventRep]
          Symmetry group representation.

      Returns
      -------
      orbit : List[OccEvent]
          The orbit of OccEvent
      )pbdoc",
      py::arg("orbit_element"), py::arg("occevent_symgroup_rep"));

  m.def(
      "make_occevent_group",
      [](occ_events::OccEvent const &occ_event,
         std::shared_ptr<occ_events::SymGroup const> const &group,
         xtal::Lattice const &lattice,
         std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep) {
        return occ_events::make_occevent_group(
            occ_event, group, lattice.lat_column_mat(), occevent_symgroup_rep);
      },
      R"pbdoc(
      Construct a subgroup which leaves an event invariant

      Parameters
      ----------
      occ_event : OccEvent
          The OccEvent that remains invariant after transformation by
          subgroup elements.

      group: List[libcasm.xtal.SymOp]
          The super group.

      lattice: xtal.Lattice
          The lattice.

      occevent_symgroup_rep: List[OccEventRep]
          Representation of `group` for transforming OccEventRep.

      Returns
      -------
      occevent_group : libcasm.sym_info.SymGroup
          The subgroup which leaves the event invariant. The head group of
          `occevent_group` is set to be the head group of `group`, which may be
          `group` itself.
      )pbdoc",
      py::arg("occ_event"), py::arg("group"), py::arg("lattice"),
      py::arg("occevent_symgroup_rep"));

  m.def(
      "make_canonical_prim_periodic_occevents",
      [](std::shared_ptr<occ_events::OccSystem const> const &system,
         clust::ClusterSpecs const &cluster_specs,
         const nlohmann::json &occevent_counter_params,
         std::vector<occ_events::OccEvent> const &custom_occevents)
          -> std::vector<occ_events::OccEvent> {
        // print errors and warnings to sys.stdout
        py::scoped_ostream_redirect redirect;
        // get canonical clusters
        auto generating_group_unitcellcoord_symgroup_rep =
            sym_info::make_unitcellcoord_symgroup_rep(
                cluster_specs.generating_group->element, *cluster_specs.prim);
        auto orbits = clust::make_prim_periodic_orbits(
            cluster_specs.prim, generating_group_unitcellcoord_symgroup_rep,
            cluster_specs.site_filter, cluster_specs.max_length,
            cluster_specs.custom_generators);
        std::vector<clust::IntegralCluster> clusters;
        for (auto const &orbit : orbits) {
          clusters.push_back(*orbit.rbegin());
        }

        // get occevent_symgroup_rep
        auto occevent_symgroup_rep = occ_events::make_occevent_symgroup_rep(
            cluster_specs.generating_group->element, *cluster_specs.prim);

        // get OccEventCounterParameters
        jsonParser json{occevent_counter_params};
        InputParser<occ_events::OccEventCounterParameters> parser(json);
        std::runtime_error error_if_invalid{
            "Error in "
            "libcasm.occ_events.make_canonical_prim_periodic_occevents"};
        report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
        occ_events::OccEventCounterParameters &params = *parser.value;

        if (json.count("print_state_info") &&
            json["print_state_info"] == true) {
          params.save_state_info = true;
          params.print_state_info =
              occ_events::OccEventCounterStateInfoPrinter(*system);
        }

        return make_prim_periodic_occevent_prototypes(
            system, clusters, occevent_symgroup_rep, params, custom_occevents);
      },
      "Documented in libcasm.occ_events._methods.py", py::arg("system"),
      py::arg("cluster_specs"), py::arg("occevent_counter_params"),
      py::arg("custom_occevents"));

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
      unitcell_index_converter : libcasm.xtal.UnitCellIndexConverter
          A :class:`~libcasm.xtal.UnitCellIndexConverter` for the supercell in which the OccEvent is located

      Returns
      -------
      (unitcell_index, equivalent_index) : tuple[int, int]
          The coordinates (unitcell_index, equivalent_index) of the
          input OccEvent. Raises if no match can be found, indicating
          the input OccEvent is not in the same orbit.
      )pbdoc",
        py::arg("occ_event"), py::arg("phenomenal_occevent"),
        py::arg("unitcell_index_converter"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
