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
  for (int i = 0; i < m_size; i++) {
    for (int j = i + 1; j < m_size; j++) {
      m_distances.push_back((cluster[i].coordinate(basicstructure) -
                             cluster[j].coordinate(basicstructure))
                                .const_cart()
                                .norm());
    }
  }
  std::sort(m_distances.begin(), m_distances.end());
}

/// \brief Construct and calculate cluster invariants,
///     including phenomenal cluster sites
ClusterInvariants::ClusterInvariants(
    IntegralCluster const &cluster, IntegralCluster const &phenomenal,
    xtal::BasicStructure const &basicstructure) {
  // save size of cluster
  m_size = cluster.size();

  // calculate distances between points
  for (int i = 0; i < m_size; i++) {
    for (int j = i + 1; j < m_size; j++) {
      m_distances.push_back((cluster[i].coordinate(basicstructure) -
                             cluster[j].coordinate(basicstructure))
                                .const_cart()
                                .norm());
    }
  }
  std::sort(m_distances.begin(), m_distances.end());

  // calculate distances between points and phenom sites
  for (int i = 0; i < cluster.size(); i++) {
    for (int j = 0; j < phenomenal.size(); j++) {
      m_phenom_distances.push_back((cluster[i].coordinate(basicstructure) -
                                    phenomenal[j].coordinate(basicstructure))
                                       .const_cart()
                                       .norm());
    }
  }
  std::sort(m_phenom_distances.begin(), m_phenom_distances.end());
}

/// \brief Number of elements in the cluster
int ClusterInvariants::size() const { return m_size; }

/// \brief const Access distances between coordinates in the cluster, sorted
/// in ascending order
const std::vector<double> &ClusterInvariants::distances() const {
  return m_distances;
}

/// \brief const Access distances between phenomenal and cluster
/// coordinates,
///     sorted in ascending order
const std::vector<double> &ClusterInvariants::phenomenal_distances() const {
  return m_phenom_distances;
}

/// \brief Check if ClusterInvariants are equal
bool almost_equal(ClusterInvariants const &A, ClusterInvariants const &B,
                  double tol) {
  return A.size() == B.size() &&
         std::equal(A.distances().cbegin(), A.distances().cend(),
                    B.distances().cbegin(),
                    [&](double a, double b) {
                      return CASM::almost_equal(a, b, tol);
                    }) &&
         std::equal(
             A.phenomenal_distances().cbegin(), A.phenomenal_distances().cend(),
             B.phenomenal_distances().cbegin(),
             [&](double a, double b) { return CASM::almost_equal(a, b, tol); });
  ;
}

/// \brief Compare ClusterInvariants
///
/// \returns True if A < B, to specified tolerance
///
/// - First compares by number of sites in cluster
/// - Then compare all distances, from longest to shortest
/// - Then compare all phenomenal - cluster distances, from longest to
/// shortest
bool compare(ClusterInvariants const &A, ClusterInvariants const &B,
             double tol) {
  // first sort by number of sites in cluster
  if (A.size() < B.size()) {
    return true;
  }
  if (A.size() > B.size()) {
    return false;
  }

  // all cluster distances
  for (int i = A.distances().size() - 1; i >= 0; i--) {
    if (CASM::almost_equal(A.distances()[i], B.distances()[i], tol)) {
      continue;
    }
    if (A.distances()[i] < B.distances()[i]) {
      return true;
    }
    if (A.distances()[i] > B.distances()[i]) {
      return false;
    }
  }

  // all cluster-phenom distances
  for (int i = A.phenomenal_distances().size() - 1; i >= 0; i--) {
    if (CASM::almost_equal(A.phenomenal_distances()[i],
                           B.phenomenal_distances()[i], tol)) {
      continue;
    }
    if (A.phenomenal_distances()[i] < B.phenomenal_distances()[i]) {
      return true;
    }
    if (A.phenomenal_distances()[i] > B.phenomenal_distances()[i]) {
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
