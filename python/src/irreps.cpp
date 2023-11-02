#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// nlohmann::json binding
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/configuration/irreps/IrrepWedge.hh"
#include "casm/configuration/irreps/VectorSpaceSymReport.hh"
#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"
#include "casm/configuration/irreps/io/json/VectorSpaceSymReport_json_io.hh"
#include "pybind11_json/pybind11_json.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_irreps, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Irreducible space decompositions

        libcasm.irreps
        --------------

        The libcasm.irreps package contains data structures and methods for finding
        irreducible space decompositions.

    )pbdoc";
  py::module::import("libcasm.sym_info");

  //
  py::class_<irreps::IrrepInfo>(m, "IrrepInfo", R"pbdoc(
            Describes an irreducible subspace.
            )pbdoc")
      .def(py::init<Eigen::MatrixXcd, Eigen::VectorXcd>(), R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          trans_mat : np.ndarray[np.float64[irrep_dim, vector_dim]]
              A shape=(`irrep_dim`,`vector_dim`) matrix that transforms a vector from
              the initial vector space into a vector in the irreducible vector space.
          characters : numpy.ndarray[numpy.complex128[m, 1]]
              A vector containing the complex character of each group operation's
              action on the irreducible vector space.
          )pbdoc",
           py::arg("trans_mat"), py::arg("characters"))
      .def_readonly("irrep_dim", &irreps::IrrepInfo::irrep_dim,
                    "int: Irreducible subspace dimension")
      .def_readonly("vector_dim", &irreps::IrrepInfo::vector_dim,
                    "int: Vector space dimension")
      .def_readonly("trans_mat", &irreps::IrrepInfo::trans_mat, R"pbdoc(
          np.ndarray[np.float64[irrep_dim, vector_dim]]: A
          shape=(`irrep_dim`,`vector_dim`) matrix that transforms a vector from the
          initial vector space into a vector in the irreducible vector space.
          )pbdoc")
      .def_readonly("characters", &irreps::IrrepInfo::trans_mat, R"pbdoc(
          numpy.ndarray[numpy.complex128[m, 1]]: A vector containing the complex
          character of each group operation's action on the irreducible vector space.
          )pbdoc")
      .def_readonly("directions", &irreps::IrrepInfo::directions,
                    R"pbdoc(
          list[list[np.ndarray[np.float64[irrep_dim,]]]: High-symmetry directions

          Vectors in the initial vector space that correspond to high-symmetry
          directions in the irreducible vector space. ``directions[i]`` is the `i`-th
          orbit of equivalent high-symmetry directions and ``len(directions[i])`` is
          the symmetric multiplicity of a direction in that orbit.
          )pbdoc")
      .def(
          "to_dict",
          [](irreps::IrrepInfo const &self) -> nlohmann::json {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the IrrepInfo as a Python dict.
          )pbdoc");

  //
  py::class_<irreps::IrrepWedge>(m, "IrrepWedge", R"pbdoc(
      A wedge in an irreducible subspace
      )pbdoc")
      .def(py::init<irreps::IrrepInfo, Eigen::MatrixXd>(), R"pbdoc(

          .. rubric: Constructor

          Parameters
          ----------
          irrep_info: IrrepInfo
              The :class:`IrrepInfo` of the irreducible subspace.
          axes: numpy.ndarray[numpy.float64[vector_dim, irrep_dim]]
              Columns of `axes` are high-symmetry directions that form the edges
              of the irreducible wedge. The number of columns equals the
              dimension of the irreducible subspace. The number of rows is usually the
              fullspace dimension, but depending on context it could be a subspace."

          )pbdoc",
           py::arg("irrep_info"), py::arg("axes"))
      .def_readonly("irrep_info", &irreps::IrrepWedge::irrep_info, R"pbdoc(
          "IrrepInfo: Irreducible space information
          )pbdoc")
      .def_readonly("axes", &irreps::IrrepWedge::axes, R"pbdoc(
          numpy.ndarray[numpy.float64[vector_dim, irrep_dim]]: High-symmetry axes

          Columns of `axes` are high-symmetry directions that form the edges
          of the irreducible wedge. The number of columns equals the
          dimension of the irreducible subspace. The number of rows is usually the
          fullspace dimension, but depending on context it could be a subspace."
          )pbdoc")
      .def_readonly(
          "mult", &irreps::IrrepWedge::mult,
          "list[int]: Symmetric multiplicity of each column of `axes`.");

  //
  py::class_<irreps::SubWedge>(m, "SubWedge", R"pbdoc(
      SubWedge is a vector of IrrepWedge, one from each irreducible subspace

      Together, the IrrepWedges that comprise the Subwedge span the entire space.
      However, it is likely that the orbit of equivalent SubWedges does not fill
      the entire space. Then multiple SubWedge combine to fill the entire space.
      )pbdoc")
      .def(py::init<std::vector<irreps::IrrepWedge> const &>(), R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          irrep_wedges: list[IrrepWedge]
              The :class:`IrrepWedge` that comprise the SubWedge
          )pbdoc",
           py::arg("irrep_wedges"))
      .def_readonly("irrep_wedges", &irreps::SubWedge::irrep_wedges, R"pbdoc(
          "list[IrrepWedge]: The :class:`IrrepWedge` that comprise the SubWedge
          )pbdoc")
      .def_readonly("trans_mat", &irreps::SubWedge::trans_mat, R"pbdoc(
          numpy.ndarray[numpy.float64[vector_dim, vector_dim]]: Transformation matrix
          to convert from a vector in terms of the SubWedge axes to a vector in the
          original vector space."
          )pbdoc");

  //
  py::class_<irreps::VectorSpaceSymReport>(m, "VectorSpaceSymReport",
                                           R"pbdoc(
      Holds results from irreducible space decompositions.

      )pbdoc")
      .def(
          py::init<std::vector<Eigen::MatrixXd>, std::vector<irreps::IrrepInfo>,
                   std::vector<irreps::SubWedge>, Eigen::MatrixXd const &,
                   std::vector<std::string>>(),
          R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          symgroup_rep: list[numpy.ndarray[numpy.float64[vector_dim, vector_dim]]]
              The symmetry group representation matrices
          irreps: list[IrrepInfo]
              The irreducible subspaces
          irreducible_wedge: list[SubWedge]
              The irreducible wedge, if calculated. Else, an empty list.
          symmetry_adapted_subspace: numpy.ndarray[numpy.float64[vector_dim, subspace_dim]]
              The symmetry adapted basis vectors, as a column vector matrix. The
              number of columns, `subspace_dim` may be less than or equal to the full
              dimension of the vector space, `vector_dim`.
          axis_glossary:
              Names given to individual axes in the initial vector space,
              corresponding to rows of `symmetry_adapted_subspace`.
          )pbdoc",
          py::arg("symgroup_rep"), py::arg("irreps"),
          py::arg("irreducible_wedge"), py::arg("symmetry_adapted_subspace"),
          py::arg("axis_glossary"))
      .def_readonly("symgroup_rep", &irreps::VectorSpaceSymReport::symgroup_rep,
                    R"pbdoc(
          list[numpy.ndarray[numpy.float64[vector_dim, vector_dim]]]: The symmetry "
          group representation matrices
          )pbdoc")
      .def_readonly("irreps", &irreps::VectorSpaceSymReport::irreps, R"pbdoc(
          list[IrrepInfo]: The irreducible subspaces
          )pbdoc")
      .def_readonly("irreducible_wedge",
                    &irreps::VectorSpaceSymReport::irreducible_wedge,
                    R"pbdoc(
          list[SubWedge]: The irreducible wedge, if calculated. Else, an empty list.
          )pbdoc")
      .def_readonly("symmetry_adapted_subspace",
                    &irreps::VectorSpaceSymReport::symmetry_adapted_subspace,
                    R"pbdoc(
          numpy.ndarray[numpy.float64[vector_dim, subspace_dim]]: The symmetry adapted
          basis vectors, as a column vector matrix. The number of columns,
          `subspace_dim` may be less than or equal to the full dimension of the vector
          space, `vector_dim`.
          )pbdoc")
      .def_readonly("axis_glossary",
                    &irreps::VectorSpaceSymReport::axis_glossary, R"pbdoc(
          list[str]: Names given to individual axes in initial vector space,
          corresponding to rows of symmetry_adapted_subspace.
          )pbdoc")
      .def_readonly("irrep_names", &irreps::VectorSpaceSymReport::irrep_names,
                    R"pbdoc(
          list[str]: Names for the irreducible subspaces

          The form of the names is `"irrep_{i}_{j}"`, where:

          - `i`: The irrep group index (starting from 1) distinguishing groups of
            identical irreps. Irreps are identical if they have the same
            character vectors.
          - `j`: The equivalent irrep index (starting from 1) distinguishing irreps
            with the same character vectors.
          )pbdoc")
      .def_readonly("irrep_axes_indices",
                    &irreps::VectorSpaceSymReport::irrep_axes_indices, R"pbdoc(
          list[list[int]]: Indices of columns in :attr:`symmetry_adapted_subspace` \
          corresponding to each irreducible subspace.

          Given a :class:`VectorSpaceSymReport`, `report`, the column index
          ``c = report.irrep_axes[i][j]`` is the column index (starting with 0)
          in ``report.symmetry_adapted_subspace`` corresponding to the `j`-th dim
          of ``report.irreps[i]`` (i.e. the `j`-th row of
          ``report.irreps[i].trans_mat``).
          )pbdoc")
      .def_readonly("irrep_wedge_axes",
                    &irreps::VectorSpaceSymReport::irrep_wedge_axes,
                    R"pbdoc(
          list[numpy.ndarray[numpy.float64[vector_dim, irrep_dim]]]: Irreducible \
          wedge axes in the full vector space, by irreducible subspace

          The matrix ``irrep_wedge_axes[i]`` are the axes for the wedge of the
          `i`-th irreducible subspace, ``irreps[i]``.
          )pbdoc")
      .def(
          "to_dict",
          [](irreps::VectorSpaceSymReport const &self) -> nlohmann::json {
            jsonParser json;
            to_json(self, json);
            return static_cast<nlohmann::json>(json);
          },
          R"pbdoc(
          Represent the VectorSpaceSymReport as a Python dict."
          )pbdoc");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
