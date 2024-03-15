#ifndef CASM_clust_IntegralCluster
#define CASM_clust_IntegralCluster

#include <set>
#include <vector>

#include "casm/configuration/clusterography/GenericCluster.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
struct UnitCellCoordRep;
}

namespace clust {

/** \defgroup Clusterography

    \brief Functions and classes related to clusters
*/

/** \defgroup IntegralCluster

    \brief Functions and classes related to IntegralCluster
    \ingroup Clusterography
    \ingroup CoordCluster
*/

/* -- IntegralCluster ------------------------------------- */

class IntegralCluster;

}  // namespace clust

/// traits, required for GenericCluster
template <>
struct traits<clust::IntegralCluster> {
  typedef xtal::UnitCellCoord Element;
  // typedef Index size_type;
  // static const std::string name;
};

namespace clust {

class IntegralCluster
    : public xtal::Translatable<GenericCluster<CRTPBase<IntegralCluster>>> {
 public:
  typedef xtal::Translatable<GenericCluster<CRTPBase<IntegralCluster>>> Base;
  using Base::Element;
  using Base::size_type;

  explicit IntegralCluster();

  explicit IntegralCluster(std::initializer_list<xtal::UnitCellCoord> elements);

  explicit IntegralCluster(std::vector<xtal::UnitCellCoord> const &elements);

  explicit IntegralCluster(std::vector<xtal::UnitCellCoord> &&elements);

  template <typename Iterator>
  IntegralCluster(Iterator begin, Iterator end);

  /// \brief Access vector of elements
  std::vector<Element> &elements();

  /// \brief const Access vector of elements
  const std::vector<Element> &elements() const;

  /// \brief Translate the cluster by a UnitCell translation
  IntegralCluster &operator+=(xtal::UnitCell trans);

 protected:
  friend GenericCluster<CRTPBase<IntegralCluster>>;

 private:
  std::vector<xtal::UnitCellCoord> m_element;
};

/// \brief Apply symmetry to IntegralCluster
IntegralCluster &apply(xtal::UnitCellCoordRep const &rep,
                       IntegralCluster &cluster);

/// \brief Apply symmetry to IntegralCluster
IntegralCluster copy_apply(xtal::UnitCellCoordRep const &rep,
                           IntegralCluster cluster);

/// \brief Convert IntegralCluster to vector of linear site
///     indices in a supercell
std::vector<Index> to_index_vector(
    IntegralCluster const &cluster,
    xtal::UnitCellCoordIndexConverter const &converter);

/// \brief Convert a set of linear site indices in a supercell
///     to an IntegralCluster
IntegralCluster cluster_from_index_vector(
    std::vector<Index> const &cluster_as_indices,
    xtal::UnitCellCoordIndexConverter const &converter);

/// \brief Convert IntegralCluster to set of linear site
///     indices in a supercell
std::set<Index> to_index_set(
    IntegralCluster const &cluster,
    xtal::UnitCellCoordIndexConverter const &converter);

/// \brief Convert a set of linear site indices in a supercell
///     to an IntegralCluster
IntegralCluster cluster_from_index_set(
    std::set<Index> const &cluster_as_indices,
    xtal::UnitCellCoordIndexConverter const &converter);

template <typename Iterator>
IntegralCluster::IntegralCluster(Iterator begin, Iterator end)
    : m_element(begin, end) {}

}  // namespace clust
}  // namespace CASM

#endif
