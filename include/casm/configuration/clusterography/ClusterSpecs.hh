#ifndef CASM_clust_ClusterSpecs
#define CASM_clust_ClusterSpecs

#include <memory>
#include <optional>
#include <set>
#include <vector>

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/IntegralClusterOrbitGenerator.hh"
#include "casm/configuration/clusterography/definitions.hh"

namespace CASM {
namespace clust {

struct ClusterSpecs {
  explicit ClusterSpecs();

  ClusterSpecs(std::shared_ptr<xtal::BasicStructure const> const &_prim,
               std::shared_ptr<SymGroup const> const &_generating_group);

  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<SymGroup const> generating_group;
  std::string site_filter_method;
  SiteFilterFunction site_filter;
  std::vector<double> max_length;
  std::vector<IntegralClusterOrbitGenerator> custom_generators;

  // --- Local clusters only ---
  std::optional<IntegralCluster> phenomenal;
  bool include_phenomenal_sites = false;
  std::vector<double> cutoff_radius;
};

// ** Filter functions **

/// \brief Generate clusters using all Site
bool all_sites_filter(const xtal::Site &site);

/// \brief Generate clusters using Site with site_occupant.size() > 1
bool alloy_sites_filter(const xtal::Site &site);

/// \brief Generate clusters using Site with specified DoF
SiteFilterFunction dof_sites_filter(const std::vector<DoFKey> &dofs = {});

/// Accept all clusters
ClusterFilterFunction all_clusters_filter();

/// Accept clusters with max pair distance less than max_length
ClusterFilterFunction max_length_cluster_filter(double max_length);

/// No sites (for null orbit, or global dof only)
CandidateSitesFunction empty_neighborhood();

/// Only sites in the origin unit cell {b, 0, 0, 0}
CandidateSitesFunction origin_neighborhood();

/// Sites within max_length distance to any site in the origin unit cell {b, 0,
/// 0, 0}
CandidateSitesFunction max_length_neighborhood(double max_length);

/// Sites within cutoff_radius distance to any site in the phenomenal cluster
CandidateSitesFunction cutoff_radius_neighborhood(
    IntegralCluster const &phenomenal, double cutoff_radius,
    bool include_phenomenal_sites = false);

}  // namespace clust
}  // namespace CASM

#endif
