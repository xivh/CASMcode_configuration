#ifndef CASM_clust_SubClusterCounter
#define CASM_clust_SubClusterCounter

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/container/Counter.hh"

namespace CASM {
namespace clust {

/// \brief Generates subclusters of a cluster
///
/// - Includes the null cluster and the original cluster
/// - Does not check uniqueness, etc.
///
class SubClusterCounter {
 public:
  /// \brief Construt with the cluster to find subclusters of
  explicit SubClusterCounter(IntegralCluster const &cluster)
      : m_cluster(cluster),
        m_site_counter(std::vector<int>(m_cluster.size(), 0),
                       std::vector<int>(m_cluster.size(), 1),
                       std::vector<int>(m_cluster.size(), 1)) {
    _set_current();
  }

  /// \brief Generate the next subcluster (if valid)
  void next() {
    ++m_site_counter;
    _set_current();
  }

  IntegralCluster const &value() const { return m_current; }

  bool valid() const { return m_site_counter.valid(); }

 private:
  void _set_current() {
    if (!valid()) {
      m_current.elements().clear();
      return;
    }
    m_current.elements().clear();
    for (Index i = 0; i < m_site_counter.size(); ++i) {
      if (m_site_counter()[i]) {
        m_current.elements().push_back(m_cluster.element(i));
      }
    }
  }

  /// the cluster we're finding subclusters of
  IntegralCluster m_cluster;

  /// The current subcluster
  IntegralCluster m_current;

  /// Indicates which sites to include (1) or not include (0) in the subcluster
  Counter<std::vector<int> > m_site_counter;
};

}  // namespace clust
}  // namespace CASM

#endif
