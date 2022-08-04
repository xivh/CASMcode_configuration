#include "casm/configuration/clusterography/IntegralCluster.hh"

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

}  // namespace clust
}  // namespace CASM
