#ifndef CASM_ClusterInvariants
#define CASM_ClusterInvariants

#include <vector>

#include "casm/configuration/clusterography/definitions.hh"

namespace CASM {
namespace clust {

class IntegralCluster;

/** \defgroup Clusterography

    \brief Functions and classes related to clusters
*/

/* -- ClusterInvariants Declarations ------------------------------------- */

/// \brief Stores cluster invariants: number of sites and site distances
///
/// Default expects (as for CoordCluster):
/// - \code ClusterType::size() \endcode
/// - \code ClusterType::coordinate(size_type index) \endcode
///
/// \ingroup Clusterography
///
class ClusterInvariants {
 public:
  /// \brief Construct and calculate cluster invariants
  ClusterInvariants(IntegralCluster const &cluster,
                    xtal::BasicStructure const &basicstructure);

  /// \brief Number of elements in the cluster
  int size() const;

  /// \brief const Access displacements between coordinates in the cluster,
  /// sorted in ascending order
  std::vector<double> const &displacement() const;

 private:
  /// \brief Number of UnitCellCoords in cluster
  int m_size;

  /// \brief Displacement between each pair of UnitCellCoords, sorted in
  /// ascending order
  std::vector<double> m_disp;
};

/// \brief Check if ClusterInvariants are equal
bool almost_equal(ClusterInvariants const &A, ClusterInvariants const &B,
                  double tol);

/// \brief Compare ClusterInvariants
bool compare(ClusterInvariants const &A, ClusterInvariants const &B,
             double tol);

/// \brief Compare std::pair<ClusterInvariants, IntegralCluster>
struct CompareCluster_f {
  typedef std::pair<ClusterInvariants, IntegralCluster> pair_type;

  CompareCluster_f(double _xtal_tol) : xtal_tol(_xtal_tol) {}

  bool operator()(pair_type const &A, pair_type const &B) const;

  double xtal_tol;
};

}  // namespace clust
}  // namespace CASM

#endif
