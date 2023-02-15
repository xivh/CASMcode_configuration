#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
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
      Data structure holding group elements and other group info, such as group-subgroup relationships.
      )pbdoc")
      .def(py::init(&make_symgroup), py::arg("element"),
           py::arg("multiplication_table"))
      .def("make_subgroup", &make_symgroup_subgroup,
           py::arg("head_group_index"), py::arg("element") = std::nullopt)
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
            std::optional<std::shared_ptr<sym_info::SymGroup const>> head_group;
            if (symgroup->head_group != nullptr) {
              head_group = symgroup->head_group;
            }
            return head_group;
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

  m.def("make_factor_group", &sym_info::make_factor_group,
        py::arg("xtal_prim"));

  m.def("make_point_group", &sym_info::make_point_group, py::arg("xtal_prim"),
        py::arg("factor_group"));

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
