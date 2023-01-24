#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/ConfigCompare.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymInfo.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// // SymGroup
//
// std::shared_ptr<config::SymGroup> make_symgroup(
//     std::vector<xtal::SymOp> const &element,
//     group::MultiplicationTable const &multiplication_table) {
//   return std::make_shared<config::SymGroup>(element, multiplication_table);
// }
//
// std::shared_ptr<config::SymGroup> make_symgroup_subgroup(
//     std::shared_ptr<config::SymGroup const> const &head_group,
//     std::set<Index> const &head_group_index,
//     std::optional<std::vector<xtal::SymOp>> element = std::nullopt) {
//   if (!element.has_value()) {
//     return std::make_shared<config::SymGroup>(head_group, head_group_index);
//   }
//   return std::make_shared<config::SymGroup>(head_group, *element,
//                                             head_group_index);
// }

// Prim

std::shared_ptr<config::Prim> make_prim(
    std::shared_ptr<xtal::BasicStructure const> const &xtal_prim) {
  return std::make_shared<config::Prim>(xtal_prim);
}

/// \brief Construct config::Prim from JSON string
std::shared_ptr<config::Prim> prim_from_json(std::string const &prim_json_str,
                                             double xtal_tol) {
  jsonParser json = jsonParser::parse(prim_json_str);
  ParsingDictionary<AnisoValTraits> const *aniso_val_dict = nullptr;
  auto basicstructure = std::make_shared<xtal::BasicStructure>(
      read_prim(json, xtal_tol, aniso_val_dict));
  return std::make_shared<config::Prim>(basicstructure);
}

/// \brief Format xtal::BasicStructure as JSON string
std::string prim_to_json(std::shared_ptr<config::Prim const> const &prim) {
  jsonParser json;
  write_prim(*prim->basicstructure, json, FRAC);
  std::stringstream ss;
  ss << json;
  return ss.str();
}

// Supercell
std::shared_ptr<config::Supercell> make_supercell(
    std::shared_ptr<config::Prim const> const &prim,
    Eigen::Matrix3l const &transformation_matrix_to_super) {
  return std::make_shared<config::Supercell>(prim,
                                             transformation_matrix_to_super);
}

// SupercellSymOp
config::SupercellSymOp make_supercell_symop(
    std::shared_ptr<config::Supercell> const &supercell,
    Index supercell_factor_group_index, std::optional<Index> translation_index,
    std::optional<Eigen::Vector3l> translation_frac,
    std::optional<Eigen::Vector3d> translation_cart) {
  if (translation_index.has_value()) {
    return config::SupercellSymOp(supercell, supercell_factor_group_index,
                                  *translation_index);
  } else if (translation_frac.has_value()) {
    xtal::UnitCell _translation_frac(*translation_frac);
    return config::SupercellSymOp(supercell, supercell_factor_group_index,
                                  _translation_frac);
  } else if (translation_cart.has_value()) {
    return config::SupercellSymOp(supercell, supercell_factor_group_index,
                                  *translation_cart);
  } else {
    throw std::runtime_error(
        "Error constructing SupercellSymOp: one of translation_index, "
        "translation_frac, or translation_cart must be provided");
  }
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_configuration, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Configuration comparison and enumeration

        libcasm.configuration
        --------------------

        The libcasm.configuration package contains data structures and methods for representing configurations, applying symmetry to configurations, and comparing configurations.

    )pbdoc";
  py::module::import("libcasm.xtal");
  py::module::import("libcasm.sym_info");

  py::class_<config::Prim, std::shared_ptr<config::Prim>>(m, "Prim", R"pbdoc(
      A data structure that includes a shared `libcasm.xtal.Prim` specifying the
      parent crystal structure and allowed degrees of freedom (DoF),
      along with the symmetry representations needed for applying symmetry to
      `libcasm.configuration.Configuration`.

      )pbdoc")
      .def(py::init(&make_prim), py::arg("xtal_prim"),
           R"pbdoc(

      .. _prim-init:

      Parameters
      ----------
      xtal_prim : libcasm.xtal.Prim
          A :class:`~libcasm.xtal.Prim`
      )pbdoc")
      .def(
          "xtal_prim",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->basicstructure;
          },
          "Returns the internal shared :class:`~libcasm.xtal.Prim`")
      .def(
          "factor_group",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->sym_info.factor_group;
          },
          "Returns the factor group.")
      .def(
          "crystal_point_group",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->sym_info.point_group;
          },
          "Returns the crystal point group.")
      .def(
          "global_dof_basis",
          [](std::shared_ptr<config::Prim const> const &prim, std::string key) {
            return prim->global_dof_info.at(key).basis();
          },
          "Returns the prim DoF basis matrix, `B`, for global DoF values of "
          "type `key`. The basis matrix converts from DoF values in the prim "
          "basis, `x_prim`, to DoF values in the standard basis, `x_standard`, "
          "according to `x_standard = B @ x_prim`.")
      .def(
          "global_dof_basis_inv",
          [](std::shared_ptr<config::Prim const> const &prim, std::string key,
             Index b) { return prim->global_dof_info.at(key).inv_basis(); },
          "Returns the prim DoF basis matrix inverse, `B_inv`, for local DoF "
          "values of type `key` on sublattice `b`. The basis matrix inverse "
          "converts from DoF values in the standard basis, `x_standard`, to "
          "DoF values in the prim basis, `x_prim`, according to `x_prim = "
          "B_inv @ x_standard`.")
      .def(
          "local_dof_basis",
          [](std::shared_ptr<config::Prim const> const &prim, std::string key,
             Index b) {
            auto const &dof_info = prim->local_dof_info.at(key);
            return dof_info[b].basis();
          },
          "Returns the prim DoF basis matrix, `B`, for local DoF values of "
          "type `key` on sublattice `b`. The basis matrix converts from DoF "
          "values in the prim basis, `x_prim`, to DoF values in the standard "
          "basis, `x_standard`, according to `x_standard = B @ x_prim`.")
      .def(
          "local_dof_basis_inv",
          [](std::shared_ptr<config::Prim const> const &prim, std::string key,
             Index b) {
            auto const &dof_info = prim->local_dof_info.at(key);
            return dof_info[b].inv_basis();
          },
          "Returns the prim DoF basis matrix inverse, `B_inv`, for local DoF "
          "values of type `key` on sublattice `b`. The basis matrix inverse "
          "converts from DoF values in the standard basis, `x_standard`, to "
          "DoF values in the prim basis, `x_prim`, according to `x_prim = "
          "B_inv @ x_standard`.")
      .def(
          "integral_site_coordinate_symgroup_rep",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->sym_info.unitcellcoord_symgroup_rep;
          },
          "Returns the symmetry group representation that describes how "
          "IntegralSiteCoordinate transform under symmetry.")
      .def(
          "occ_symgroup_rep",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->sym_info.occ_symgroup_rep;
          },
          "Returns the symmetry group representation that describes how "
          "occupant indices transform under symmetry.")
      .def(
          "atom_position_symgroup_rep",
          [](std::shared_ptr<config::Prim const> const &prim) {
            return prim->sym_info.atom_position_symgroup_rep;
          },
          "Returns the symmetry group representation that describes how atom "
          "position indices transform under symmetry.")
      .def_static(
          "from_json", &prim_from_json,
          "Construct a Prim from a JSON-formatted string. The `Prim reference "
          "<https://prisms-center.github.io/CASMcode_docs/formats/casm/"
          "crystallography/BasicStructure/>`_ documents the expected JSON "
          "format.",
          py::arg("prim_json_str"), py::arg("xtal_tol") = TOL)
      .def("to_json", &prim_to_json,
           "Represent the Prim as a JSON-formatted string. The `Prim reference "
           "<https://prisms-center.github.io/CASMcode_docs/formats/casm/"
           "crystallography/BasicStructure/>`_ documents the expected JSON "
           "format.");

  py::class_<config::Supercell, std::shared_ptr<config::Supercell>>(
      m, "Supercell", R"pbdoc(
      A data structure that stores the supercell transformation matrix and
      the symmetry representations needed for applying symmetry to
      `libcasm.configuration.Configuration` associated with a particular
      supercell.

      )pbdoc")
      .def(py::init(&make_supercell), py::arg("prim"),
           py::arg("transformation_matrix_to_super"),
           R"pbdoc(Constructor

      Parameters
      ----------
      prim : libcasm.configuration.Prim
          A :class:`~libcasm.configuration.Prim`
      transformation_matrix_to_super : array_like, shape=(3,3), dtype=int
          The transformation matrix, T, relating the superstructure lattice vectors, S, to the unit structure lattice vectors, L, according to S = L @ T, where S and L are shape=(3,3)  matrices with lattice vectors as columns.
      )pbdoc")
      .def(
          "prim",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->prim;
          },
          "Returns the internal shared :class:`~libcasm.configuration.Prim`")
      .def(
          "xtal_prim",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->prim->basicstructure;
          },
          "Returns the internal shared :class:`~libcasm.xtal.Prim`")
      .def(
          "transformation_matrix_to_super",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->superlattice.transformation_matrix_to_super();
          },
          "Returns the supercell transformation matrix.")
      .def(
          "superlattice",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->superlattice.superlattice();
          },
          "Returns the supercell lattice.")
      .def(
          "prim_lattice",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->superlattice.prim_lattice();
          },
          "Returns the prim lattice.")
      .def(
          "factor_group",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->sym_info.factor_group;
          },
          "Returns the supercell factor group, which is the subgroup of the "
          "prim factor group that leaves the supercell lattice vectors "
          "invariant.")
      .def(
          "factor_group_permutations",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->sym_info.factor_group_permutations;
          },
          "Returns the factor group permutations, where "
          "`factor_group_permutations()[i]` describes how "
          "`factor_group().element(i)` permutes supercell site DoF values. "
          "When permuting site occupants, the following convention is used, "
          "`after[l] = before[permutation[l]]`. When applying symmetry to "
          "anisotropic DoF values, the convention used by CASM, (for example, "
          "by `apply_supercell_sym_op`), is that the DoF values are "
          "transformed first, then permuted.")
      .def(
          "translation_permutations",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->sym_info.translation_permutations;
          },
          "Returns the translation permutations, where "
          "`translations_permutations()[i]` describes how the translation "
          "vector corresponding to translating the origin unit cell to the "
          "unitcell specified by linear index `l`, permutes supercell site DoF "
          "values. When permuting site occupants, the following convention is "
          "used, `after[l] = before[permutation[l]]`.")
      .def(
          "n_sites",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->unitcellcoord_index_converter.total_sites();
          },
          "Returns the number of sites in the supercell.")
      .def(
          "n_unitcells",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            return supercell->unitcell_index_converter.total_sites();
          },
          "Returns the unit cells in the supercell.")
      .def(
          "n_occupants",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            std::vector<Index> n_occupants;
            auto const &converter = supercell->unitcellcoord_index_converter;
            auto const &xtal_prim = supercell->prim->basicstructure;
            auto const &unique_names = xtal_prim->unique_names();
            Index n_sites = converter.total_sites();
            for (Index l = 0; l < n_sites; ++l) {
              Index b = converter(l).sublattice();
              n_occupants.push_back(unique_names[b].size());
            }
            return n_occupants;
          },
          "Returns a `List[int]`, such that `n_occupants[l]` is the number of "
          "occupants allowed on site `l`.")
      .def(
          "occ_dof",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            std::vector<std::vector<std::string>> occ_dof;
            auto const &converter = supercell->unitcellcoord_index_converter;
            auto const &xtal_prim = supercell->prim->basicstructure;
            auto const &unique_names = xtal_prim->unique_names();
            Index n_sites = converter.total_sites();
            for (Index l = 0; l < n_sites; ++l) {
              Index b = converter(l).sublattice();
              occ_dof.push_back(unique_names[b]);
            }
            return occ_dof;
          },
          "Returns a `List[List[str]]`, such that `occ_dof[l][s]` is the name "
          "('orientation name') of the occupant corresponding to occupation "
          "value `s` on site `l`.")
      .def(
          "coordinate_cart",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            auto const &xtal_prim = supercell->prim->basicstructure;
            Index n_sites = converter.total_sites();
            Eigen::MatrixXd R(3, n_sites);
            for (Index l = 0; l < n_sites; ++l) {
              xtal::UnitCellCoord bijk = converter(l);
              R.col(l) = bijk.coordinate(*xtal_prim).const_cart();
            }
            return R;
          },
          "Returns the basis site positions, as columns of a matrix, in "
          "Cartesian coordinates")
      .def(
          "coordinate_frac",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            auto const &xtal_prim = supercell->prim->basicstructure;
            Index n_sites = converter.total_sites();
            Eigen::MatrixXd R(3, n_sites);
            for (Index l = 0; l < n_sites; ++l) {
              xtal::UnitCellCoord bijk = converter(l);
              R.col(l) = bijk.coordinate(*xtal_prim).const_frac();
            }
            return R;
          },
          "Returns the basis site positions, as columns of a matrix, in "
          "fractional coordinates of the prim lattice vectors.")
      .def(
          "sublattice_indices",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            std::vector<Index> sublattice_indices;
            auto const &converter = supercell->unitcellcoord_index_converter;
            Index n_sites = converter.total_sites();
            for (Index l = 0; l < n_sites; ++l) {
              sublattice_indices.push_back(converter(l).sublattice());
            }
            return sublattice_indices;
          },
          "Returns the sublattice indices, as a List[int], of each site "
          "in the supercell.")
      .def(
          "unitcell_indices",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            Index n_sites = converter.total_sites();
            Eigen::MatrixXl R(3, n_sites);
            for (Index l = 0; l < n_sites; ++l) {
              xtal::UnitCellCoord bijk = converter(l);
              R.col(l) = bijk.unitcell();
            }
            return R;
          },
          "Returns the integer unit cell indices, as columns of a matrix, of "
          "each site in the supercell.")
      .def(
          "linear_unitcell_indices",
          [](std::shared_ptr<config::Supercell const> const &supercell) {
            std::vector<Index> unitcell_indices;
            auto const &xtal_prim = supercell->prim->basicstructure;
            Index N_sublat = xtal_prim->basis().size();
            Index N_unitcells =
                supercell->unitcell_index_converter.total_sites();
            for (Index b = 0; b < N_sublat; ++b) {
              for (Index i = 0; i < N_unitcells; ++i) {
                unitcell_indices.push_back(i);
              }
            }
            return unitcell_indices;
          },
          "Returns the linear unitcell index for each site in the supercell.")
      .def(
          "linear_unitcell_index",
          [](std::shared_ptr<config::Supercell const> const &supercell,
             Eigen::Vector3l const &unitcell) {
            return supercell->unitcell_index_converter(unitcell);
          },
          py::arg("unitcell"),
          "Returns the linear unitcell index for the given unit cell "
          "coordinates.")
      .def(
          "linear_site_index",
          [](std::shared_ptr<config::Supercell const> const &supercell,
             xtal::UnitCellCoord const &integral_site_coordinate) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            return converter(integral_site_coordinate);
          },
          py::arg("integral_site_coordinate"),
          "Converts integral_site_coordinate to the linear site index of the "
          "periodic equivalent site within the supercell.")
      .def(
          "within",
          [](std::shared_ptr<config::Supercell const> const &supercell,
             xtal::UnitCellCoord const &integral_site_coordinate) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            return converter.bring_within(integral_site_coordinate);
          },
          py::arg("integral_site_coordinate"),
          "Returns the periodic equivalent integral_site_coordinate within the "
          "supercell.")
      .def(
          "integral_site_coordinate",
          [](std::shared_ptr<config::Supercell const> const &supercell,
             Index l) {
            auto const &converter = supercell->unitcellcoord_index_converter;
            return converter(l);
          },
          "Returns the integral_site_coordinate of a site in the supercell.")
      .def(py::self < py::self,
           "Sorts supercells by size then how canonical the lattice vectors "
           "are. Only supercells with the same prim can be compared.")
      .def(py::self <= py::self,
           "Sorts supercells by size then how canonical the lattice vectors "
           "are. Only supercells with the same prim can be compared.")
      .def(py::self > py::self,
           "Sorts supercells by size then how canonical the lattice vectors "
           "are. Only supercells with the same prim can be compared.")
      .def(py::self >= py::self,
           "Sorts supercells by size then how canonical the lattice vectors "
           "are. Only supercells with the same prim can be compared.")
      .def(py::self == py::self,
           "True if supercells are equal. Only supercells with the same prim "
           "can be compared.")
      .def(py::self != py::self,
           ""
           "True if supercells are not equal. Only supercells with the same "
           "prim can be compared.");

  m.def(
      "is_canonical_supercell",
      [](std::shared_ptr<config::Supercell const> const &supercell) {
        return is_canonical(*supercell);
      },
      py::arg("init_supercell"),
      "Return true if supercell lattice is in canonical form");

  m.def(
      "make_canonical_supercell",
      [](std::shared_ptr<config::Supercell const> const &supercell) {
        return make_canonical_form(*supercell);
      },
      py::arg("init_supercell"),
      "Return a supercell that compares greater to all equivalent supercells "
      "generated using the prim crystal point group");

  m.def(
      "to_canonical_supercell",
      [](std::shared_ptr<config::Supercell const> const &supercell) {
        return to_canonical(*supercell);
      },
      py::arg("supercell"),
      "Return the :class:`~libcasm.xtal.SymOp` that makes a supercell lattice "
      "canonical");

  m.def(
      "from_canonical_supercell",
      [](std::shared_ptr<config::Supercell const> const &supercell) {
        return from_canonical(*supercell);
      },
      py::arg("supercell"),
      "Return the :class:`~libcasm.xtal.SymOp` that makes a supercell lattice "
      "from the canonical supercell lattice");

  m.def(
      "make_equivalent_supercells",
      [](std::shared_ptr<config::Supercell const> const &supercell) {
        return make_equivalents(*supercell);
      },
      py::arg("supercell"),
      "Return a list of supercells with distinct, but symmetrically equivalent "
      "lattice vectors, using the prim crystal point group.");

  // SupercellSymOp -- declare class
  py::class_<config::SupercellSymOp> pySupercellSymOp(m, "SupercellSymOp",
                                                      R"pbdoc(
      Represents and allows iteration over symmetry operations consistent
      with a given Supercell, combining pure factor group and pure
      translation operations.

    )pbdoc");

  // Configuration -- declare class
  py::class_<config::Configuration> pyConfiguration(m, "Configuration", R"pbdoc(
      A data structure encapsulating configuration DoF values and all
      information necessary to apply supercell factor group operations.

    )pbdoc");

  // SupercellSymOp -- define functions
  pySupercellSymOp
      .def(py::init(&make_supercell_symop), py::arg("supercell"),
           py::arg("supercell_factor_group_index"),
           py::arg("translation_index"), py::arg("translation_frac"),
           py::arg("translation_cart"),
           R"pbdoc(
      Construct a representation that allows applying symmetry to
      configurations in a particular supercell.

      Notes
      -----
      Only one of `translation_index`, `translation_frac`, or
      `translation_cart` should be provided.

      Parameters
      ----------
      supercell : libcasm.configuration.Supercell
          The supercell defining the periodicity of the configuration to
          be transformed.
      supercell_factor_group_index : int
          Index into the supercell factor group.
      translation_index : int, optional
          Linear index of translation to apply. Equivalent to linear
          unitcell_index.
      translation_frac : array_like, shape=(3,), dtype=int, optional
          Translation, as multiples of prim lattice vectors.
      translation_cart : array_like, shape=(3,), dtype=double, optional
          Translation, in Cartesian coordinates.
      )pbdoc")
      .def("supercell", &config::SupercellSymOp::supercell,
           "Return the supercell.")
      .def(
          "next", [](config::SupercellSymOp &op) { ++op; },
          "Transform this to represent the next operation in the supercell "
          "factor group.")
      .def_static("begin", &config::SupercellSymOp::begin,
                  "Construct a representation for the first operation in the "
                  "supercell factor group.")
      .def_static("end", &config::SupercellSymOp::end,
                  "Construct a representation for the past-the-last operation "
                  "in the supercell factor group.")
      .def_static("translation_begin",
                  &config::SupercellSymOp::translation_begin,
                  "Construct a representation for the first operation in the "
                  "set of supercell translations.")
      .def_static("translation_end", &config::SupercellSymOp::translation_end,
                  "Construct a representation for the past-the-last operation "
                  "in the set of supercell translations.")
      .def("supercell_factor_group_index",
           &config::SupercellSymOp::supercell_factor_group_index,
           "Index into the supercell factor group of the element applied by "
           "this supercell factor group operation.")
      .def("prim_factor_group_index",
           &config::SupercellSymOp::prim_factor_group_index,
           "Index into the prim factor group of the element applied by this "
           "supercell factor group operation.")
      .def("translation_index", &config::SupercellSymOp::translation_index,
           "Index of the translation applied by this supercell factor group "
           "operation.")
      .def("translation_frac", &config::SupercellSymOp::translation_frac,
           "The translation applied by this supercell factor group operation, "
           "in fractional coordinates relative to the prim lattice vectors.")
      .def("permute_index", &config::SupercellSymOp::permute_index,
           py::arg("i"),
           "Returns the index of the site containing the site DoF values that "
           "will be permuted onto site i")
      .def("to_symop", &config::SupercellSymOp::to_symop,
           "Return the SymOp for the current operation")
      .def("combined_permute", &config::SupercellSymOp::combined_permute,
           "Returns the combination of factor group operation permutation and "
           "translation permutation")
      .def("inverse", &config::SupercellSymOp::inverse,
           "Returns the inverse supercell operation")
      .def(
          "__mul__",
          [](config::SupercellSymOp const &self,
             config::SupercellSymOp const &rhs) { return self * rhs; },
          py::arg("rhs"),
          "Return the supercell operation equivalent to applying first rhs and "
          "then this.")
      .def(py::self < py::self,
           "Sorts SupercellSymOp. Only SupercellSymOp with the same supercell "
           "can be compared.")
      .def(py::self <= py::self,
           "Sorts SupercellSymOp. Only SupercellSymOp with the same supercell "
           "can be compared.")
      .def(py::self > py::self,
           "Sorts SupercellSymOp. Only SupercellSymOp with the same supercell "
           "can be compared.")
      .def(py::self >= py::self,
           "Sorts SupercellSymOp. Only SupercellSymOp with the same supercell "
           "can be compared.")
      .def(py::self == py::self,
           "True if SupercellSymOp are equal. Only SupercellSymOp with the "
           "same supercell can be compared.")
      .def(py::self != py::self,
           "True if SupercellSymOp are not equal. Only SupercellSymOp with the "
           "same supercell can be compared.")
      .def(
          "__mul__",
          [](config::SupercellSymOp const &self,
             config::Configuration const &configuration) {
            return copy_apply(self, configuration);
          },
          py::arg("configuration"),
          "Creates a copy of the configuration and applies the symmetry "
          "operation represented by this SupercellSymOp to its degree of "
          "freedom "
          "(DoF) values.")
      .def(
          "__mul__",
          [](config::SupercellSymOp const &self,
             xtal::UnitCellCoord const &integral_site_coordinate) {
            return copy_apply(self, integral_site_coordinate);
          },
          py::arg("integral_site_coordinate"),
          "Creates a copy of `integral_site_coordinate` and applies the "
          "symmetry operation represented by this SupercellSymOp.");

  // Configuration -- define functions
  pyConfiguration
      .def(py::init<std::shared_ptr<config::Supercell const> const &>(),
           py::arg("supercell"),
           R"pbdoc(
      Construct a default-valued Configuration

      Parameters
      ----------
      supercell : libcasm.configuration.Supercell
          The supercell defining the periodicity of the configuration.
      )pbdoc")
      .def(
          "supercell",
          [](config::Configuration const &configuration) {
            return configuration.supercell;
          },
          "Returns the internal shared "
          ":class:`~libcasm.configuration.Supercell`")
      .def(
          "occupation",
          [](config::Configuration const &configuration) {
            return configuration.dof_values.occupation;
          },
          "Returns the site occupation values, as indices into the allowed "
          "occupants on the corresponding basis site.")
      .def(
          "set_occupation",
          [](config::Configuration &configuration,
             Eigen::VectorXi const &occupation) {
            if (occupation.size() !=
                configuration.dof_values.occupation.size()) {
              throw std::runtime_error(
                  "Error in set_occupation: occupation vector size may not be "
                  "changed");
            }
            return configuration.dof_values.occupation = occupation;
          },
          "Sets the site occupation values. Changing the size of the "
          "occupation vector results in an exception.")
      .def(
          "occ",
          [](config::Configuration const &configuration, Index l) {
            return configuration.dof_values.occupation(l);
          },
          py::arg("l"),
          "Returns the site occupation value on the specified site.")
      .def(
          "set_occ",
          [](config::Configuration &configuration, Index l, int s) {
            configuration.dof_values.occupation(l) = s;
          },
          py::arg("l"), py::arg("s"),
          "Set the site occupation value on site `l` to value `s`. Values are "
          "not checked for validity, but `l` is expected to be in range `[0, "
          "supercell.n_sites())`, `s` is expected to be in range `[0, "
          "supercell.n_occupants()[l])`.")
      .def(
          "global_dof_values",
          [](config::Configuration const &configuration, std::string key) {
            return configuration.dof_values.global_dof_values.at(key);
          },
          py::arg("key"),
          "Returns global DoF values of type `key`, in the prim basis.")
      .def(
          "global_standard_dof_values",
          [](config::Configuration const &configuration, std::string key) {
            auto const &prim = configuration.supercell->prim;
            return clexulator::global_to_standard_values(
                configuration.dof_values.global_dof_values.at(key),
                prim->global_dof_info.at(key));
          },
          py::arg("key"),
          "Returns global DoF values of type `key`, in the standard basis.")
      .def(
          "set_global_dof_values",
          [](config::Configuration &configuration, std::string key,
             Eigen::VectorXd const &dof_values) {
            Eigen::VectorXd &_dof_values =
                configuration.dof_values.global_dof_values.at(key);
            if (_dof_values.size() != dof_values.size()) {
              throw std::runtime_error(
                  "Error in set_global_dof_values: size may not be changed");
            }
            return _dof_values = dof_values;
          },
          py::arg("key"), py::arg("dof_values"),
          "Set global DoF values of type `key`, in the prim basis.")
      .def(
          "set_global_standard_dof_values",
          [](config::Configuration &configuration, std::string key,
             Eigen::VectorXd const &standard_dof_values) {
            auto const &prim = configuration.supercell->prim;
            Eigen::VectorXd &_dof_values =
                configuration.dof_values.global_dof_values.at(key);
            _dof_values = clexulator::global_from_standard_values(
                standard_dof_values, prim->global_dof_info.at(key));
          },
          py::arg("key"), py::arg("standard_dof_values"),
          "Set global DoF values of type `key`, in the standard basis.")
      .def(
          "local_dof_values",
          [](config::Configuration const &configuration, std::string key) {
            return configuration.dof_values.local_dof_values.at(key);
          },
          py::arg("key"),
          "Returns local DoF values of type `key`, in the prim basis.")
      .def(
          "local_standard_dof_values",
          [](config::Configuration const &configuration, std::string key) {
            auto const &supercell = configuration.supercell;
            auto const &prim = supercell->prim;
            auto const &xtal_prim = prim->basicstructure;
            Index N_sublat = xtal_prim->basis().size();
            Index N_unitcells =
                supercell->unitcell_index_converter.total_sites();
            return clexulator::local_to_standard_values(
                configuration.dof_values.local_dof_values.at(key), N_sublat,
                N_unitcells, prim->local_dof_info.at(key));
          },
          py::arg("key"),
          "Returns local DoF values of type `key`, in the standard basis.")
      .def(
          "set_local_dof_values",
          [](config::Configuration &configuration, std::string key,
             Eigen::MatrixXd const &dof_values) {
            Eigen::MatrixXd &_dof_values =
                configuration.dof_values.local_dof_values.at(key);
            if (_dof_values.rows() != dof_values.rows()) {
              throw std::runtime_error(
                  "Error in set_local_dof_values: number of rows may not be "
                  "changed");
            }
            if (_dof_values.cols() != dof_values.cols()) {
              throw std::runtime_error(
                  "Error in set_local_dof_values: number of cols may not be "
                  "changed");
            }
            return _dof_values = dof_values;
          },
          py::arg("key"), py::arg("dof_values"),
          "Set local DoF values of type `key`, in the prim basis.")
      .def(
          "set_local_standard_dof_values",
          [](config::Configuration &configuration, std::string key,
             Eigen::MatrixXd const &standard_dof_values) {
            auto const &supercell = configuration.supercell;
            auto const &prim = supercell->prim;
            auto const &xtal_prim = prim->basicstructure;
            Index N_sublat = xtal_prim->basis().size();
            Index N_unitcells =
                supercell->unitcell_index_converter.total_sites();
            Eigen::MatrixXd &_dof_values =
                configuration.dof_values.local_dof_values.at(key);
            _dof_values = clexulator::local_from_standard_values(
                standard_dof_values, N_sublat, N_unitcells,
                prim->local_dof_info.at(key));
          },
          py::arg("key"), py::arg("standard_dof_values"),
          "Set local DoF values of type `key`, in the standard basis.")
      .def(
          "local_dof_site_value",
          [](config::Configuration const &configuration, std::string key,
             Index l) {
            return configuration.dof_values.local_dof_values.at(key).col(l);
          },
          py::arg("key"), py::arg("l"),
          "Returns local DoF values of type `key`, in the prim basis, on site "
          "`l`.")
      .def(
          "set_local_dof_site_value",
          [](config::Configuration &configuration, std::string key, Index l,
             Eigen::VectorXd const &site_dof_value) {
            configuration.dof_values.local_dof_values.at(key).col(l) =
                site_dof_value;
          },
          py::arg("key"), py::arg("l"), py::arg("site_dof_value"),
          "Set the local DoF values of type `key`, in the prim basis, on site "
          "`l`.")
      .def(
          "local_standard_dof_site_value",
          [](config::Configuration const &configuration, std::string key,
             Index l) {
            auto const &prim = configuration.supercell->prim;
            auto const &converter =
                configuration.supercell->unitcellcoord_index_converter;
            Index b = converter(l).sublattice();
            auto const &dof_info = prim->local_dof_info.at(key);
            auto const &M = configuration.dof_values.local_dof_values.at(key);
            return dof_info[b].basis() * M.col(l).head(dof_info[b].dim());
          },
          py::arg("key"), py::arg("l"),
          "Returns local DoF values of type `key`, in the standard basis, on "
          "site `l`.")
      .def(
          "set_local_standard_dof_site_value",
          [](config::Configuration &configuration, std::string key, Index l,
             Eigen::VectorXd const &standard_site_dof_value) {
            auto const &prim = configuration.supercell->prim;
            auto const &converter =
                configuration.supercell->unitcellcoord_index_converter;
            Index b = converter(l).sublattice();
            auto const &dof_info = prim->local_dof_info.at(key);
            auto &M = configuration.dof_values.local_dof_values.at(key);
            M.col(l).head(dof_info[b].dim()) =
                dof_info[b].inv_basis() * standard_site_dof_value;
          },
          py::arg("key"), py::arg("l"), py::arg("standard_site_dof_value"),
          "Set the local DoF values of type `key`, in the standard basis, on "
          "site `l`.")
      .def(py::self < py::self,
           "Sorts configurations, first by supercell, then global DoF, then "
           "occupation DoF, then local continuous DoF. Only configurations "
           "with the same prim can be compared.")
      .def(py::self <= py::self,
           "Sorts configurations, first by supercell, then global DoF, then "
           "occupation DoF, then local continuous DoF. Only configurations "
           "with the same prim can be compared.")
      .def(py::self > py::self,
           "Sorts configurations, first by supercell, then global DoF, then "
           "occupation DoF, then local continuous DoF. Only configurations "
           "with the same prim can be compared.")
      .def(py::self >= py::self,
           "Sorts configurations, first by supercell, then global DoF, then "
           "occupation DoF, then local continuous DoF. Only configurations "
           "with the same prim can be compared.")
      .def(py::self == py::self,
           "True if configurations are equal, or approximately equal up the "
           "lattice tolerance if there continuous DoF. Only configurations "
           "with the same prim can be compared.")
      .def(py::self != py::self,
           "True if configurations are equal, or approximately equal up the "
           "lattice tolerance if there continuous DoF. Only configurations "
           "with the same prim can be compared.");

  m.def(
      "apply",
      [](config::SupercellSymOp const &op,
         config::Configuration &configuration) {
        return apply(op, configuration);
      },
      py::arg("supercell_symop"), py::arg("configuration"),
      "Applies the symmetry operation represented by the SupercellSymOp to "
      "transform the configuration's degree of freedom (DoF) values.");

  m.def(
      "copy_apply",
      [](config::SupercellSymOp const &op,
         config::Configuration const &configuration) {
        return copy_apply(op, configuration);
      },
      py::arg("supercell_symop"), py::arg("configuration"),
      "Creates a copy of the configuration and applies the symmetry "
      "operation represented by this SupercellSymOp to its degree of freedom "
      "(DoF) values.");

  m.def(
      "apply",
      [](config::SupercellSymOp const &op,
         xtal::UnitCellCoord &integral_site_coordinate) {
        return apply(op, integral_site_coordinate);
      },
      py::arg("supercell_symop"), py::arg("integral_site_coordinate"),
      "Applies the symmetry operation represented by the SupercellSymOp to "
      "transform `integral_site_coordinate`.");

  m.def(
      "copy_apply",
      [](config::SupercellSymOp const &op,
         xtal::UnitCellCoord const &integral_site_coordinate) {
        return copy_apply(op, integral_site_coordinate);
      },
      py::arg("supercell_symop"), py::arg("configuration"),
      "Creates a copy of `integral_site_coordinate` and applies the symmetry "
      "operation represented by this SupercellSymOp");

  m.def(
      "is_canonical_configuration",
      [](config::Configuration const &configuration,
         std::optional<config::SupercellSymOp> begin,
         std::optional<config::SupercellSymOp> end) {
        auto const &supercell = configuration.supercell;
        if (!begin.has_value()) {
          begin = config::SupercellSymOp::begin(supercell);
        }
        if (!end.has_value()) {
          end = config::SupercellSymOp::end(supercell);
        }
        return is_canonical(configuration, *begin, *end);
      },
      py::arg("configuration"), py::arg("begin") = std::nullopt,
      py::arg("end") = std::nullopt,
      "Return true if configuration is in canonical form, given the provided "
      "SupercellSymOp (default is the supercell factor group).");

  m.def(
      "make_canonical_configuration",
      [](config::Configuration const &configuration,
         bool in_canonical_supercell,
         std::optional<std::vector<config::SupercellSymOp>> subgroup) {
        if (in_canonical_supercell) {
          if (subgroup.has_value()) {
            throw std::runtime_error(
                "Error in make_canonical_configuration: "
                "`in_canonical_supercell` may not be used in combination with "
                "`subgroup`");
          }
          return make_in_canonical_supercell(configuration);
        } else if (subgroup.has_value()) {
          return make_canonical_form(configuration, subgroup->begin(),
                                     subgroup->end());
        } else {
          auto const &supercell = configuration.supercell;
          auto begin = config::SupercellSymOp::begin(supercell);
          auto end = config::SupercellSymOp::end(supercell);
          return make_canonical_form(configuration, begin, end);
        }
      },
      py::arg("configuration"), py::arg("in_canonical_supercell") = false,
      py::arg("subgroup") = std::nullopt,
      R"pbdoc(
      Return the canonical configuration

      The canonical configuration is the configuration which compares greater
      than all equivalent supercells given a symmetry group. By default, the
      supercell factor group is used. The parameter
      `in_canonical_supercell` can be used to specify that the supercell
      should also be made canonical. A subgroup of the supercell factor
      group may be used through the `subgroup` parameter, but not in combination
      with `in_canonical_supercell`.

      Parameters
      ----------
      configuration : libcasm.configuration.Configuration
          The initial configuration.
      in_canonical_supercell : bool, default=False
          If True, the canonical configuration is found in the
          canonical supercell. Otherwise the supercell is not changed.
      subgroup : List[libcasm.configuration.SupercellSymOp], optional
          If provided, the canonical configuration will be found with
          respect to a subgroup of the supercell factor group instead of
          the complete supercell factor group. These operations are supercell
          specifici, so this parameter may not be used with
          `in_canonical_supercell == True`.
      )pbdoc");

  m.def(
      "to_canonical_configuration",
      [](config::Configuration const &configuration,
         std::optional<std::vector<config::SupercellSymOp>> subgroup) {
        if (subgroup.has_value()) {
          return to_canonical(configuration, subgroup->begin(),
                              subgroup->end());
        } else {
          auto const &supercell = configuration.supercell;
          auto begin = config::SupercellSymOp::begin(supercell);
          auto end = config::SupercellSymOp::end(supercell);
          return to_canonical(configuration, begin, end);
        }
      },
      py::arg("configuration"), py::arg("subgroup") = std::nullopt,
      "Return an operation that makes a configuration canonical with respect "
      "to the supercell factor group (default) or a subgroup of the supercell "
      "factor group.");

  m.def(
      "from_canonical_configuration",
      [](config::Configuration const &configuration,
         std::optional<std::vector<config::SupercellSymOp>> subgroup) {
        if (subgroup.has_value()) {
          return from_canonical(configuration, subgroup->begin(),
                                subgroup->end());
        } else {
          auto const &supercell = configuration.supercell;
          auto begin = config::SupercellSymOp::begin(supercell);
          auto end = config::SupercellSymOp::end(supercell);
          return from_canonical(configuration, begin, end);
        }
      },
      py::arg("configuration"), py::arg("subgroup") = std::nullopt,
      "Return an operation that makes a configuration from the canonical "
      "configuration with respect to the supercell factor group (default) or a "
      "subgroup of the supercell factor group.");

  m.def(
      "make_invariant_subgroup",
      [](config::Configuration const &configuration,
         std::optional<std::set<Index>> site_indices,
         std::optional<std::vector<config::SupercellSymOp>> group) {
        if (!site_indices.has_value()) {
          if (group.has_value()) {
            return make_invariant_subgroup(configuration, group->begin(),
                                           group->end());
          } else {
            auto const &supercell = configuration.supercell;
            auto begin = config::SupercellSymOp::begin(supercell);
            auto end = config::SupercellSymOp::end(supercell);
            return make_invariant_subgroup(configuration, begin, end);
          }
        } else {
          if (group.has_value()) {
            return make_invariant_subgroup(configuration, site_indices.value(),
                                           group->begin(), group->end());
          } else {
            auto const &supercell = configuration.supercell;
            auto begin = config::SupercellSymOp::begin(supercell);
            auto end = config::SupercellSymOp::end(supercell);
            return make_invariant_subgroup(configuration, site_indices.value(),
                                           begin, end);
          }
        }
      },
      py::arg("configuration"), py::arg("site_indices") = std::nullopt,
      py::arg("group") = std::nullopt,
      "Return the subgroup (as a List[libcasm.configuration.SupercellSymOp]) "
      "that leaves a configuration invariant. If `site_indices` are provided "
      "(a set of linear index of sites in the supercell), the subgroup is "
      "restricted such that it does not mix the given sites and other sites. "
      "By default, the subgroup is found with respect to the supercell factor "
      "group. Optionally, another `group` (itself a subgroup of the supercell "
      "factor group) may be provided.");

  m.def(
      "make_invariant_subgroup",
      [](std::set<Index> site_indices,
         std::optional<std::vector<config::SupercellSymOp>> group) {
        if (group.has_value()) {
          return make_invariant_subgroup(site_indices, group->begin(),
                                         group->end());
        } else {
          auto supercell = group->begin()->supercell();
          auto begin = config::SupercellSymOp::begin(supercell);
          auto end = config::SupercellSymOp::end(supercell);
          return make_invariant_subgroup(site_indices, begin, end);
        }
      },
      py::arg("site_indices"), py::arg("group") = std::nullopt,
      "Return the subgroup (as a List[libcasm.configuration.SupercellSymOp]) "
      "that does not mix the given sites (a set of linear index of sites in "
      "the supercell) and other sites. By default, the subgroup is found with "
      "respect to the supercell factor group. Optionally, another `group` "
      "(itself a subgroup of the supercell factor group) may be provided.");

  m.def(
      "make_equivalent_configurations",
      [](config::Configuration const &configuration,
         std::optional<std::vector<config::SupercellSymOp>> subgroup) {
        if (subgroup.has_value()) {
          return make_equivalents(configuration, subgroup->begin(),
                                  subgroup->end());
        } else {
          auto const &supercell = configuration.supercell;
          auto begin = config::SupercellSymOp::begin(supercell);
          auto end = config::SupercellSymOp::end(supercell);
          return make_equivalents(configuration, begin, end);
        }
      },
      py::arg("configuration"), py::arg("subgroup") = std::nullopt,
      "Return the distinct symmetrically equivalent configurations, with "
      "respect to the supercell factor group (default) or a subgroup of the "
      "supercell factor group.");

  m.def(
      "copy_configuration",
      [](config::Configuration const &motif,
         std::shared_ptr<config::Supercell const> const &supercell,
         Eigen::Vector3l origin) {
        return copy_configuration(motif, supercell, xtal::UnitCell(origin));
      },
      py::arg("motif"), py::arg("supercell"),
      py::arg("origin") = Eigen::Vector3l::Zero(),
      R"pbdoc(
      Copy motif configuration DoF values into a supercell

      Parameters
      ----------
      motif : libcasm.configuration.Configuration
          The initial configuration, with DoF values to be filled into the supercell.
      supercell : libcasm.configuration.Supercell
          The supercell to be filled by the motif configuration.
      origin : array_like of int, shape=(3,)
          The UnitCell indicating which unit cell in the initial configuration
          is the origin in new configuration
      )pbdoc");

  m.def(
      "copy_transformed_configuration",
      [](Index prim_factor_group_index, Eigen::Vector3l const &translation,
         config::Configuration const &motif,
         std::shared_ptr<config::Supercell const> const &supercell,
         Eigen::Vector3l origin) {
        return copy_configuration(prim_factor_group_index, translation, motif,
                                  supercell, xtal::UnitCell(origin));
      },
      py::arg("prim_factor_group_index"), py::arg("translation"),
      py::arg("motif"), py::arg("supercell"),
      py::arg("origin") = Eigen::Vector3l::Zero(),
      R"pbdoc(
      Copy transformed motif configuration DoF values into a supercell

      Copies DoF values as if `motif` is transformed by the prim factor group
      operation with index `factor_group_index`, then translated by
      `translation`, then copied starting from `origin`. In other words, sites
      map according to:

          new_config_unitcellcoord + origin = fg * motif_unitcellcoord + trans


      Parameters
      ----------
      prim_factor_group_index : int
          Index of prim factor group operation which transforms the initial
          configuration
      translation : array_like of int, shape=(3,)
          Lattice translation applied after the prim factor group operation.
      motif : libcasm.configuration.Configuration
          The initial configuration, with DoF values to be filled into the supercell.
      supercell : libcasm.configuration.Supercell
          The supercell to be filled by the motif configuration.
      origin : array_like of int, shape=(3,)
          The UnitCell indicating which unit cell in the initial configuration
          is the origin in new configuration
      )pbdoc");

  m.def("make_all_super_configurations", &config::make_all_super_configurations,
        py::arg("motif"), py::arg("supercell"),
        R"pbdoc(
      Make all equivalent configurations with respect to the prim factor group
      that fill a supercell

      Parameters
      ----------
      motif : libcasm.configuration.Configuration
          The initial configuration, with DoF values to be filled into the supercell.
      supercell : libcasm.configuration.Supercell
          The supercell to be filled by the motif configuration.

      Returns
      -------
      super_configurations : List[libcasm.configuration.Configuration]
          All configurations equivalent with respect to the prim factor
          group which fit in the given supercell.
      )pbdoc");

  m.def("is_primitive_configuration", &config::is_primitive,
        py::arg("configuration"),
        "Return true if no translations within the supercell result in the "
        "same configuration");

  m.def("make_primitive_configuration", &config::make_primitive,
        py::arg("configuration"),
        "Return the primitive configuration. Does not apply any symmetry "
        "operations. Use `make_canonical_configuration` with "
        "`in_canonical_supercell=True` aftwards to obtain the "
        "primitive canonical configuration in the canonical supercell.");

  m.def("make_global_dof_matrix_rep", &config::make_global_dof_matrix_rep,
        py::arg("group"), py::arg("key"),
        "Make the matrix representation of `group` that describes the "
        "transformation of a particular global DoF, specified by `key`");

  m.def("make_local_dof_matrix_rep", &config::make_local_dof_matrix_rep,
        py::arg("group"), py::arg("key"), py::arg("site_indices"),
        "Make the matrix representation of `group` that describes the "
        "transformation of a particular local DoF, specified by `key` "
        "(may be \"occ\") amongst a subset of supercell sites, "
        "specified by `site_indices` (a set of linear index of sites in "
        "the supercell).");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
