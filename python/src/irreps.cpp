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
#include "pybind11_json/pybind11_json.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

irreps::VectorSpaceSymReport make_VectorSpaceSymReport(
    std::vector<Eigen::MatrixXd> const &symgroup_rep,
    std::vector<irreps::IrrepInfo> const &irreps,
    std::vector<irreps::SubWedge> const &irreducible_wedge,
    Eigen::MatrixXd const &symmetry_adapted_subspace,
    std::vector<std::string> const &axis_glossary) {
  irreps::VectorSpaceSymReport report;
  report.symgroup_rep = symgroup_rep;
  report.irreps = irreps;
  report.irreducible_wedge = irreducible_wedge;
  report.symmetry_adapted_subspace = symmetry_adapted_subspace;
  report.axis_glossary = axis_glossary;
  return report;
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_irreps, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Irreducible space decompositions

        libcasm.irreps
        --------------

        The libcasm.irreps package contains data structures and methods for finding irreducible space decompositions.

    )pbdoc";
  py::module::import("libcasm.sym_info");

  //
  py::class_<irreps::IrrepInfo>(m, "IrrepInfo", R"pbdoc(
            Describes an irreducible subspace.
            )pbdoc")
      .def(py::init<Eigen::MatrixXcd, Eigen::VectorXcd>(), R"pbdoc(
          Construct IrrepInfo.

          trans_mat : np.ndarray
              An irrep_dim() x vector_dim() matrix that transforms a vector from the
              initial vector space into a vector in the irreducible vector space.
          characters : np.ndarray
              A vector containing the complex character of each group operation's
              action on the irreducible vector space.
          )pbdoc",
           py::arg("trans_mat"), py::arg("characters"))
      .def_readonly("irrep_dim", &irreps::IrrepInfo::irrep_dim,
                    "int: Irreducible subspace dimension")
      .def_readonly("vector_dim", &irreps::IrrepInfo::vector_dim,
                    "int: Vector space dimension")
      .def_readonly(
          "trans_mat", &irreps::IrrepInfo::trans_mat,
          "np.ndarray: An irrep_dim() x vector_dim() matrix that transforms "
          "a vector from the initial vector space into a vector in the "
          "irreducible vector space")
      .def_readonly(
          "characters", &irreps::IrrepInfo::trans_mat,
          "np.ndarray: An irrep_dim() x vector_dim() matrix that transforms "
          "a vector from the initial vector space into a vector in the "
          "irreducible vector space.")
      .def_readonly("directions", &irreps::IrrepInfo::directions,
                    R"pbdoc(
      list[list[np.ndarray]]: High-symmetry directions
          Vectors in the initial vector space that correspond to high-symmetry
          directions in the irreducible vector space. ``directions[i]`` is the `i`-th
          orbit of equivalent high-symmetry directions and ``len(directions[i])`` is
          the symmetric multiplicity of a direction in that orbit.
      )pbdoc");

  //
  py::class_<irreps::IrrepWedge>(m, "IrrepWedge", R"pbdoc(
      A wedge in an irreducible subspace
      )pbdoc")
      .def(py::init<irreps::IrrepInfo, Eigen::MatrixXd>(), R"pbdoc(
          Construct IrrepWedge
          )pbdoc",
           py::arg("irrep_info"), py::arg("axes"))
      .def_readonly(
          "irrep_info", &irreps::IrrepWedge::irrep_info,
          ":class:`~libcasm.irreps.IrrepInfo`: Irreducible space information")
      .def_readonly("axes", &irreps::IrrepWedge::axes, R"pbdoc(
          np.ndarray: High-symmetry axes
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
          Construct SubWedge
          )pbdoc",
           py::arg("irrep_wedges"))
      .def_readonly("irrep_wedges", &irreps::SubWedge::irrep_wedges,
                    "list[:class:`~libcasm.irreps.IrrepWedge`]: The IrrepWedge "
                    "that comprise the SubWedge")
      .def_readonly("trans_mat", &irreps::SubWedge::trans_mat, R"pbdoc(
          np.ndarray:
              Transformation matrix to convert from a vector in terms of the SubWedge
              axes to a vector in the original vector space."
          )pbdoc");

  //
  py::class_<irreps::VectorSpaceSymReport>(m, "VectorSpaceSymReport",
                                           R"pbdoc(
      Holds results from :func:`~libcasm.irreps.vector_space_sym_report`.

      )pbdoc")
      .def(py::init(&make_VectorSpaceSymReport), R"pbdoc(
           Construct a VectorSpaceSymReport.
           )pbdoc",
           py::arg("symgroup_rep"), py::arg("irreps"),
           py::arg("irreducible_wedge"), py::arg("symmetry_adapted_subspace"),
           py::arg("axis_glossary"))
      .def_readonly(
          "symgroup_rep", &irreps::VectorSpaceSymReport::symgroup_rep,
          "list[np.ndarray]: The symmetry group representation matrices")
      .def_readonly(
          "irreps", &irreps::VectorSpaceSymReport::irreps,
          "list[:class:`~libcasm.irreps.IrrepInfo`]: The irreducible subspaces")
      .def_readonly(
          "irreducible_wedge", &irreps::VectorSpaceSymReport::irreducible_wedge,
          "list[:class:`~libcasm.irreps.SubWedge`]: The irreducible wedge "
          "(if calculated)")
      .def_readonly("symmetry_adapted_subspace",
                    &irreps::VectorSpaceSymReport::symmetry_adapted_subspace,
                    "np.ndarray: The symmetry adapted basis vectors, as a "
                    "column vector matrix")
      .def_readonly(
          "axis_glossary", &irreps::VectorSpaceSymReport::axis_glossary,
          "list[str]: Names given to individual axes in initial (un-adapted) "
          "vector "
          "space, corresponding to rows of symmetry_adapted_subspace.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
