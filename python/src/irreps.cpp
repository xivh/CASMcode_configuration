#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// nlohmann::json binding
#define JSON_USE_IMPLICIT_CONVERSIONS 0
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/subgroups.hh"
#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/configuration/irreps/IrrepWedge.hh"
#include "casm/configuration/irreps/VectorSpaceSymReport.hh"
#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"
#include "casm/configuration/irreps/io/json/VectorSpaceSymReport_json_io.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "pybind11_json/pybind11_json.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

namespace CASM {
namespace irreps {
typedef group::Group<Eigen::MatrixXd> MatrixRepGroup;
}
}  // namespace CASM

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

std::shared_ptr<irreps::MatrixRepGroup> make_matrixrepgroup(
    std::vector<Eigen::MatrixXd> const &elements, double abs_tol) {
  return std::make_shared<irreps::MatrixRepGroup>(group::make_group(
      elements,
      [](Eigen::MatrixXd const &A,
         Eigen::MatrixXd const &B) -> Eigen::MatrixXd { return A * B; },
      [=](Eigen::MatrixXd const &A, Eigen::MatrixXd const &B) -> bool {
        return CASM::almost_equal(A, B, abs_tol);
      }));
}

std::shared_ptr<irreps::MatrixRepGroup const> make_matrixrepgroup_subgroup(
    std::shared_ptr<irreps::MatrixRepGroup const> const &group,
    std::set<Index> const &head_group_index) {
  std::shared_ptr<irreps::MatrixRepGroup const> head_group;
  if (!group->head_group) {
    head_group = group;
  } else {
    head_group = group->head_group;
  }
  return std::make_shared<irreps::MatrixRepGroup>(head_group, head_group_index);
}

irreps::IrrepDecomposition make_IrrepDecomposition(
    irreps::MatrixRep const &matrix_rep,
    std::optional<irreps::GroupIndices> head_group,
    std::optional<Eigen::MatrixXd> init_subspace, bool allow_complex,
    double abs_tol) {
  if (matrix_rep.size() == 0) {
    throw std::runtime_error(
        "Error in make_IrrepDecomposition: matrix_rep.size() == 0");
  }
  for (auto const &M : matrix_rep) {
    if (M.rows() != matrix_rep[0].rows()) {
      throw std::runtime_error(
          "Error in make_IrrepDecomposition: matrix reps are different sizes");
    }
    if (M.rows() != M.cols()) {
      throw std::runtime_error(
          "Error in make_IrrepDecomposition: matrix_rep is not square");
    }
  }
  if (!head_group.has_value()) {
    irreps::GroupIndices _head_group;
    Index i = 0;
    for (auto const &rep : matrix_rep) {
      _head_group.insert(i);
      ++i;
    }
    head_group = _head_group;
  }
  if (!init_subspace.has_value()) {
    int dim = matrix_rep[0].rows();
    init_subspace = Eigen::MatrixXd::Identity(dim, dim);
  }

  std::shared_ptr<irreps::MatrixRepGroup> matrixrepgroup =
      make_matrixrepgroup(matrix_rep, abs_tol);

  std::function<irreps::GroupIndicesOrbitSet()> make_cyclic_subgroups_f =
      [=]() { return group::make_cyclic_subgroups(*matrixrepgroup); };
  std::function<irreps::GroupIndicesOrbitSet()> make_all_subgroups_f = [=]() {
    return group::make_all_subgroups(*matrixrepgroup);
  };

  std::optional<Log> log;
  return irreps::IrrepDecomposition(matrix_rep, *head_group, *init_subspace,
                                    make_cyclic_subgroups_f,
                                    make_all_subgroups_f, allow_complex, log);
}

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
  py::class_<irreps::MatrixRepGroup, std::shared_ptr<irreps::MatrixRepGroup>>(
      m, "MatrixRepGroup", R"pbdoc(
      Data structure holding group elements and other group info, such as
      group-subgroup relationships.

      .. rubric:: Special methods

      - Get the `i`-th element of the group using: ``element = group[i]``

      )pbdoc")
      .def(py::init(&make_matrixrepgroup), R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------

          elements: list[numpy.ndarray[numpy.float64[m, m]]
              The matrix representation of elements of the group.
          abs_tol: float = :data:`~libcasm.casmglobal.TOL`
              The absolute tolerance, used to construct the group multiplication
              table.
          )pbdoc",
           py::arg("elements"), py::arg("abs_tol") = CASM::TOL)
      .def("make_subgroup", &make_matrixrepgroup_subgroup, R"pbdoc(
          Make a subgroup

          Parameters
          ----------
          head_group_index: list[int]
              Indices of elements in the head group (which may or may not be
              `self`) of elements to include in subgroup.
          )pbdoc",
           py::arg("head_group_index"))
      .def_property_readonly(
          "elements",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) { return matrixrepgroup->element; },
          R"pbdoc(
          list[numpy.ndarray[numpy.float64[m, m]]]: The elements list
          )pbdoc")
      .def(
          "__getitem__",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup,
             Index i) { return matrixrepgroup->element[i]; },
          py::arg("i"),
          R"pbdoc(
          numpy.ndarray[numpy.float64[m, m]]: The `i`-th element in the group.
          )pbdoc")
      .def_property_readonly(
          "multiplication_table",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) {
            return matrixrepgroup->multiplication_table;
          },
          R"pbdoc(
          list[list[int]]: The multiplication table.

          The multiplication table element `multiplication_table[i][j] == k`
          represents that
          ``np.allclose(elements[k], elements[i] @ elements[j]) == True``.
          )pbdoc")
      .def(
          "conjugacy_classes",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) {
            return make_conjugacy_classes(*matrixrepgroup);
          },
          R"pbdoc(
          Returns the conjugacy classes

          Returns
          -------
          conjugacy_classes: list[list[int]]
              ``conjugacy_classes[i]`` is a list of the indices of elements in
              the `i`-th class."
          )pbdoc")
      .def(
          "mult",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup,
             Index i, Index j) { return matrixrepgroup->mult(i, j); },
          py::arg("i"), py::arg("j"),
          R"pbdoc(
          Returns the index of the element product.

          Parameters
          ----------
          i: int
              lhs element index.
          j: int
              rhs element index.

          Returns
          -------
          k: int
              The index ``k == multiplication_table[i][j]``.
          )pbdoc")
      .def_property_readonly(
          "is_subgroup",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) {
            return matrixrepgroup->head_group != nullptr;
          },
          R"pbdoc(
          bool: True if this is a subgroup, False otherwise.
          )pbdoc")
      .def_property_readonly(
          "head_group",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) {
            std::optional<std::shared_ptr<irreps::MatrixRepGroup const>>
                head_group;
            if (matrixrepgroup->head_group != nullptr) {
              head_group = matrixrepgroup->head_group;
            }
            return head_group;
          },
          R"pbdoc(
          Optional[MatrixRepGroup]: The head group if this is a subgroup, else None.
          )pbdoc")
      .def_property_readonly(
          "head_group_index",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) { return matrixrepgroup->head_group_index; },
          R"pbdoc(
          list[int]: The list of head group indices (guaranteed sorted)

          If this is the head group, then
          ``group.head_group_index == [0, 1, 2, ...]``.

          If this is a sub group, then ``subgroup.element[i]`` is the same element
          as ``subgroup.head_group.element[subgroup.head_group_index[i]]``.

          )pbdoc")
      .def_property_readonly(
          "inverse_index",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup) { return matrixrepgroup->inverse_index; },
          R"pbdoc(
          list[int]: The list of inverse indices

          Represents that ``group.element[group.inverse_index[i]]`` is the
          inverse of ``group.element[i]``.
          )pbdoc")
      .def(
          "inv",
          [](std::shared_ptr<irreps::MatrixRepGroup const> const
                 &matrixrepgroup,
             Index i) { return matrixrepgroup->inv(i); },
          py::arg("i"), "Return the index of the inverse element.",
          R"pbdoc(
          Returns the index of the inverse of an element

          Parameters
          ----------
          i: int
              The element index.

          Returns
          -------
          i_inverse: int
              The index the inverse of the `i`-th element.
          )pbdoc");

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
          np.ndarray[np.complex128[irrep_dim, vector_dim]]: A
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
          matrix_rep: list[numpy.ndarray[numpy.float64[vector_dim, vector_dim]]]
              The group element representation matrices
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
          py::arg("matrix_rep"), py::arg("irreps"),
          py::arg("irreducible_wedge"), py::arg("symmetry_adapted_subspace"),
          py::arg("axis_glossary"))
      .def_readonly("matrix_rep", &irreps::VectorSpaceSymReport::symgroup_rep,
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

  //
  py::class_<irreps::IrrepDecomposition>(m, "IrrepDecomposition", R"pbdoc(
            Performs irreducible subspace construction and symmetrization
            )pbdoc")
      .def(py::init<>(&make_IrrepDecomposition), R"pbdoc(

          .. rubric:: Constructor

          Parameters
          ----------
          matrix_rep: list[np.ndarray[np.float64]]
              Full space matrix representation
          head_group: Optional[set[int]] = None
              Group used to find irreps, as indices into `matrix_rep`. If
              None, the entire group is used.
          init_subspace: Optional[np.ndarray[np.float64]] = None
              Space in which to find irreducible subspaces. This space is
              expanded, if necessary, by application of `matrix_rep` and
              orthogonalization to form an invariant subspace (i.e. column
              space does not change upon application of elements in head_group).
              If None, use the identity matrix of dimension matching
              `matrix_rep`.
          allow_complex: bool = True
              If True, all irreps may be complex-valued, if False, complex
              irreps are combined to form real representations
          abs_tol: float = :data:`~libcasm.casmglobal.TOL`
              The absolute tolerance, used to construct a group multiplication
              table.
          )pbdoc",
           py::arg("matrix_rep"), py::arg("head_group") = std::nullopt,
           py::arg("init_subspace") = std::nullopt,
           py::arg("allow_complex") = true, py::arg("abs_tol") = CASM::TOL)
      .def_readonly("matrix_rep", &irreps::IrrepDecomposition::fullspace_rep,
                    "Full space matrix representation")
      .def_readonly("head_group", &irreps::IrrepDecomposition::head_group,
                    "Group used to find irreps, as indices into `matrix_rep`")
      .def_readonly("subspace", &irreps::IrrepDecomposition::subspace,
                    "Subspace, invariant with respect to `head_group` in which "
                    "irreps are found.")
      .def_readonly("symmetry_adapted_subspace",
                    &irreps::IrrepDecomposition::symmetry_adapted_subspace,
                    "Symmetry adapted subspace, formed from `irreps`.")
      .def_readonly("irreps", &irreps::IrrepDecomposition::irreps, R"pbdoc(
          Irreducible spaces, symmetrized to align the irreducible space bases
          along high symmetry directions. Irreps are found in the `subspace`
          and then converted to full space dimension (meaning
          `irrep[i].vector_dim() == full space dimension` and
          `sum_i irrep[i].irrep_dim() == subspace columns`).
          )pbdoc")
      .def(
          "make_symmetry_report",
          [](irreps::IrrepDecomposition const &self, bool calc_wedges,
             std::optional<std::vector<std::string>> glossary)
              -> irreps::VectorSpaceSymReport {
            if (!glossary.has_value()) {
              Index dim = self.fullspace_rep[0].rows();
              std::vector<std::string> _glossary;
              for (Index i = 0; i < dim; ++i) {
                _glossary.push_back(std::string("x") + std::to_string(i + 1));
              }
              glossary = _glossary;
            }

            return irreps::vector_space_sym_report(self, calc_wedges,
                                                   *glossary);
          },
          R"pbdoc(
          Construct a VectorSpaceSymReport

          Parameters
          ----------
          calc_wedges : bool = False
              If True, calculate the irreducible wedges for the vector space.
              This may take a long time, but provides the symmetrically unique
              portions of the vector space, which is useful for enumeration.
          glossary: Optional[list[str]] = None
              If provided, a description of each dimension of the vector space.
          )pbdoc",
          py::arg("calc_wedges") = false, py::arg("glossary") = std::nullopt);

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
