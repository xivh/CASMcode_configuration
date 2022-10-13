#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/SymType.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

// SymGroup

std::shared_ptr<sym_info::SymGroup> make_symgroup(
    std::vector<xtal::SymOp> const &element,
    group::MultiplicationTable const &multiplication_table) {
  return std::make_shared<sym_info::SymGroup>(element, multiplication_table);
}

std::shared_ptr<sym_info::SymGroup> make_symgroup_subgroup(
    std::shared_ptr<sym_info::SymGroup const> const &head_group,
    std::set<Index> const &head_group_index,
    std::optional<std::vector<xtal::SymOp>> element = std::nullopt) {
  if (!element.has_value()) {
    return std::make_shared<sym_info::SymGroup>(head_group, head_group_index);
  }
  return std::make_shared<sym_info::SymGroup>(head_group, *element,
                                              head_group_index);
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_sym_info, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Crystallographic symmetry information

        libcasm.sym_info
        ----------------

        The libcasm.sym_info package constructs crystallographic symmetry information.

    )pbdoc";
  py::module::import("libcasm.xtal");

  py::class_<sym_info::SymGroup, std::shared_ptr<sym_info::SymGroup>>(
      m, "SymGroup",
      R"pbdoc(
      Holds group elements and group info.
      )pbdoc")
      .def(py::init(&make_symgroup), py::arg("element"),
           py::arg("multiplication_table"),
           R"pbdoc(

      Parameters
      ----------
      element : List[libcasm.xtal.SymOp]
          List of elements of the symmetry group.
      multiplication_table : List[List[int]]
          The group multiplication table, where `k == multiplication_table[i][j]`
          indicates that element[k] == element[i] * element[j]. This must be
          provide for a head group, but it is determined for a subgroup.
      )pbdoc")
      .def("make_subgroup", &make_symgroup_subgroup,
           py::arg("head_group_index"), py::arg("element") = std::nullopt,
           R"pbdoc(
      Construct a subgroup

      Parameters
      ----------
      head_group_index : Set[int]
          Indices of elements to include in the subgroup.
      element : List[libcasm.xtal.SymOp], optional
          List of elements in the subgroup. This is optional, and if not provided
          the head group elements corresponding to the head_group_index will be
          used. It may be provided to specify particular elements for a subgroup
          of a factor group, such as the elements of a cluster group. Must be
          listed in ascending head group index order.
      )pbdoc")
      .def(
          "elements",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->element;
          },
          "Returns the elements list.")
      .def(
          "__getitem__",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup,
             Index i) { return symgroup->element[i]; },
          py::arg("i"), "Returns an element from group.")
      .def(
          "multiplication_table",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->multiplication_table;
          },
          "Returns the multiplication table.")
      .def(
          "conjugacy_classes",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return make_conjugacy_classes(*symgroup);
          },
          "Returns a List[List[int]], where `conjugacy_classes[i]` is a list "
          "of the indices of elements in class 'i'")
      .def(
          "mult",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup, Index i,
             Index j) { return symgroup->mult(i, j); },
          py::arg("i"), py::arg("j"),
          "Returns the index `k == multiplication_table[i][j]`.")
      .def(
          "is_subgroup",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->head_group != nullptr;
          },
          "Returns True if this is a subgroup.")
      .def(
          "head_group",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->head_group;
          },
          "If this is a subgroup, return the head group, else return None.")
      .def(
          "head_group_index",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->head_group_index;
          },
          "Return the list of head group indices")
      .def(
          "inverse_index",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup) {
            return symgroup->inverse_index;
          },
          "Return the list of head group indices")
      .def(
          "inv",
          [](std::shared_ptr<sym_info::SymGroup const> const &symgroup,
             Index i) { return symgroup->inv(i); },
          py::arg("i"), "Return the index of the inverse element.");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
