#include "casm/configuration/clusterography/ClusterInvariants.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace clust {

/// \brief Construct and calculate cluster invariants
ClusterInvariants::ClusterInvariants(
    IntegralCluster const &cluster,
    xtal::BasicStructure const &basicstructure) {
  // save size of cluster
  m_size = cluster.size();

  // calculate distances between points
  std::vector<double> disp;
  for (int i = 0; i < m_size; i++) {
    for (int j = i + 1; j < m_size; j++) {
      m_disp.push_back((cluster[i].coordinate(basicstructure) -
                        cluster[j].coordinate(basicstructure))
                           .const_cart()
                           .norm());
    }
  }
  std::sort(m_disp.begin(), m_disp.end());
}

/// \brief Number of elements in the cluster
int ClusterInvariants::size() const { return m_size; }

/// \brief const Access displacements between coordinates in the cluster, sorted
/// in ascending order
const std::vector<double> &ClusterInvariants::displacement() const {
  return m_disp;
}

/// \brief Check if ClusterInvariants are equal
bool almost_equal(ClusterInvariants const &A, ClusterInvariants const &B,
                  double tol) {
  return A.size() == B.size() &&
         std::equal(A.displacement().cbegin(), A.displacement().cend(),
                    B.displacement().cbegin(), [&](double a, double b) {
                      return CASM::almost_equal(a, b, tol);
                    });
}

/// \brief Compare ClusterInvariants
///
/// \returns True if A < B, to specified tolerance
///
/// - First compares by number of sites in cluster
/// - Then compare all displacements, from longest to shortest
///
bool compare(ClusterInvariants const &A, ClusterInvariants const &B,
             double tol) {
  // first sort by number of sites in cluster
  if (A.size() < B.size()) {
    return true;
  }
  if (A.size() > B.size()) {
    return false;
  }

  // all displacements
  for (int i = A.displacement().size() - 1; i >= 0; i--) {
    if (CASM::almost_equal(A.displacement()[i], B.displacement()[i], tol)) {
      continue;
    }
    if (A.displacement()[i] < B.displacement()[i]) {
      return true;
    }
    if (A.displacement()[i] > B.displacement()[i]) {
      return false;
    }
  }
  return false;
}

bool CompareCluster_f::operator()(pair_type const &A,
                                  pair_type const &B) const {
  if (compare(A.first, B.first, xtal_tol)) {
    return true;
  }
  if (compare(B.first, A.first, xtal_tol)) {
    return false;
  }
  return A.second < B.second;
}

}  // namespace clust
}  // namespace CASM
