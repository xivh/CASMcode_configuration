#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

/// CASM - Python binding code
namespace CASMpy {

using namespace CASM;

xtal::UnitCellCoord site_from_list(std::vector<Index> const &site) {
  if (site.size() != 4) {
    throw std::runtime_error(
        "Error integral sit coordinate from list: list size != 4");
  }
  return xtal::UnitCellCoord(site[0], site[1], site[2], site[3]);
}

std::vector<Index> site_to_list(xtal::UnitCellCoord const &site) {
  std::vector<Index> list;
  for (Index i = 0; i < 4; ++i) {
    list.push_back(site[i]);
  }
  return list;
}

clust::IntegralCluster &append_site(clust::IntegralCluster &cluster,
                                    xtal::UnitCellCoord const &site) {
  cluster.elements().push_back(site);
  return cluster;
}

clust::IntegralCluster &append_site_list(clust::IntegralCluster &cluster,
                                         std::vector<Index> const &site) {
  cluster.elements().push_back(site_from_list(site));
  return cluster;
}

clust::IntegralCluster make_cluster_from_list(
    std::vector<std::vector<Index>> const &list) {
  clust::IntegralCluster cluster;
  for (auto const &site : list) {
    append_site_list(cluster, site);
  }
  return cluster;
}

}  // namespace CASMpy

PYBIND11_DECLARE_HOLDER_TYPE(T, std::shared_ptr<T>);

PYBIND11_MODULE(_clusterography, m) {
  using namespace CASMpy;

  m.doc() = R"pbdoc(
        Cluster classes and methods

        libcasm.clusterography
        --------------------

        The libcasm.clusterography package contains data structures and methods for generating clusters, orbits of clusters, and neighborhoods.

    )pbdoc";
  py::module::import("libcasm.xtal");

  py::class_<clust::IntegralCluster>(m, "Cluster",
                                     R"pbdoc(
      A cluster of IntegralSiteCoordinate
      )pbdoc")
      .def(py::init<std::vector<xtal::UnitCellCoord> const &>(),
           py::arg("sites"),
           R"pbdoc(

      Parameters
      ----------
      sites : List[libcasm.xtal.IntegralSiteCoordinate]
          List of sites in the cluster. May be empty.
      )pbdoc")
      .def_static("from_list", &make_cluster_from_list, py::arg("sites"),
                  R"pbdoc(

      Parameters
      ----------
      sites : List[List[int]]
          List of sites (as [b, i, j, k]) in the cluster. May be empty.
      )pbdoc")
      .def(
          "sites",
          [](clust::IntegralCluster const &cluster) {
            return cluster.elements();
          },
          "Returns the cluster sites")
      .def("append", &append_site, "Append to cluster sites")
      .def("append", &append_site_list, "Append to cluster sites")
      .def(
          "to_list",
          [](clust::IntegralCluster const &cluster) {
            std::vector<std::vector<Index>> list;
            for (auto const &site : cluster.elements()) {
              list.push_back(site_to_list(site));
            }
            return list;
          },
          "Returns the cluster sites as a list of sites (list of [b, i, j, k])")
      .def(
          "site",
          [](clust::IntegralCluster const &cluster, Index i) {
            return cluster.element(i);
          },
          "Returns the i-th site in the cluster")
      .def(
          "size",
          [](clust::IntegralCluster const &cluster) {
            return cluster.elements().size();
          },
          "Returns number of sites in the cluster")
      .def(
          "__getitem__",
          [](clust::IntegralCluster const &cluster, Index i) {
            return cluster.element(i);
          },
          "Returns the i-th site in the cluster")
      .def(
          "__add__",
          [](clust::IntegralCluster const &cluster,
             Eigen::Vector3l const &translation) {
            return cluster + translation;
          },
          "Translate a cluster by adding unit cell indices")
      .def(
          "__iadd__",
          [](clust::IntegralCluster &cluster,
             Eigen::Vector3l const &translation) {
            return cluster += translation;
          },
          "Translate a cluster by adding unit cell indices")
      .def(
          "__sub__",
          [](clust::IntegralCluster const &cluster,
             Eigen::Vector3l const &translation) {
            return cluster - translation;
          },
          "Translate a cluster by subtracting unit cell indices")
      .def(
          "__isub__",
          [](clust::IntegralCluster &cluster,
             Eigen::Vector3l const &translation) {
            return cluster -= translation;
          },
          "Translate a cluster by subtracting unit cell indices");

#ifdef VERSION_INFO
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
