#include "casm/configuration/clusterography/IntegralCluster.hh"

#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace clust {

IntegralCluster::IntegralCluster() {}

IntegralCluster::IntegralCluster(
    std::initializer_list<xtal::UnitCellCoord> elements)
    : m_element(elements) {}

IntegralCluster::IntegralCluster(
    std::vector<xtal::UnitCellCoord> const &elements)
    : m_element(elements) {}

IntegralCluster::IntegralCluster(std::vector<xtal::UnitCellCoord> &&elements)
    : m_element(std::move(elements)) {}

/// \brief Access vector of elements
std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() {
  return m_element;
}

/// \brief const Access vector of elements
const std::vector<xtal::UnitCellCoord> &IntegralCluster::elements() const {
  return m_element;
}

/// \brief Translate the cluster by a UnitCell translation
IntegralCluster &IntegralCluster::operator+=(xtal::UnitCell trans) {
  for (auto it = this->begin(); it != this->end(); ++it) {
    *it += trans;
  }
  return *this;
}

/// \brief Apply symmetry to IntegralCluster
IntegralCluster &apply(xtal::UnitCellCoordRep const &rep,
                       IntegralCluster &cluster) {
  for (auto &unitcellcoord : cluster) {
    apply(rep, unitcellcoord);
  }
  return cluster;
}

/// \brief Apply symmetry to IntegralCluster
IntegralCluster copy_apply(xtal::UnitCellCoordRep const &rep,
                           IntegralCluster cluster) {
  apply(rep, cluster);
  return cluster;
}

/// \brief Convert IntegralCluster to vector of linear site
///     indices in a supercell
std::vector<Index> to_index_vector(
    IntegralCluster const &cluster,
    xtal::UnitCellCoordIndexConverter const &converter) {
  std::vector<Index> index_vector;
  for (auto const &site : cluster) {
    index_vector.push_back(converter(site));
  }
  return index_vector;
}

/// \brief Convert a set of linear site indices in a supercell
///     to an IntegralCluster
IntegralCluster cluster_from_index_vector(
    std::vector<Index> const &cluster_as_indices,
    xtal::UnitCellCoordIndexConverter const &converter) {
  IntegralCluster cluster;
  for (auto linear_site_index : cluster_as_indices) {
    cluster.elements().push_back(converter(linear_site_index));
  }
  return cluster;
}

/// \brief Convert IntegralCluster to set of linear site
///     indices in a supercell
std::set<Index> to_index_set(
    IntegralCluster const &cluster,
    xtal::UnitCellCoordIndexConverter const &converter) {
  std::set<Index> index_set;
  for (auto const &site : cluster) {
    index_set.insert(converter(site));
  }
  return index_set;
}

/// \brief Convert a set of linear site indices in a supercell
///     to an IntegralCluster
IntegralCluster cluster_from_index_set(
    std::set<Index> const &cluster_as_indices,
    xtal::UnitCellCoordIndexConverter const &converter) {
  IntegralCluster cluster;
  for (auto linear_site_index : cluster_as_indices) {
    cluster.elements().push_back(converter(linear_site_index));
  }
  return cluster;
}

}  // namespace clust
}  // namespace CASM
