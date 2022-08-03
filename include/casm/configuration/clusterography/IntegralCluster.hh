#ifndef CASM_clust_IntegralCluster
#define CASM_clust_IntegralCluster

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
                           IntegralCluster &cluster);

template <typename Iterator>
IntegralCluster::IntegralCluster(Iterator begin, Iterator end)
    : m_element(begin, end) {}

}  // namespace clust
}  // namespace CASM

#endif
