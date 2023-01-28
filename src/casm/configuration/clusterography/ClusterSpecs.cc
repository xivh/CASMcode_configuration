#include "casm/configuration/clusterography/ClusterSpecs.hh"

#include "casm/configuration/clusterography/ClusterInvariants.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/container/Counter.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace clust {

namespace {

/// \brief Output the neighborhood of UnitCellCoord within max_radius of any
/// site in unit cell
///
/// \param unit The unit cell xtal::BasicStructure
/// \param max_radius The neighborhood distance cutoff
/// \param site_filter A filter function that returns true for CoordType that
///        should be considered for the neighborhood
/// \param result Output iterator for container of UnitCellCoord
/// \param xtal_tol Crystallography tolerance used to contstruct UnitCellCoord
/// from CoordType
///
/// \returns Output iterator after generating the neighborhood
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator>
OutputIterator neighborhood(xtal::BasicStructure const &unit, double max_radius,
                            SiteFilterFunction site_filter,
                            OutputIterator result, double xtal_tol) {
  auto dim = unit.lattice().enclose_sphere(max_radius);
  EigenCounter<Eigen::Vector3i> grid_count(-dim, dim,
                                           Eigen::Vector3i::Constant(1));
  xtal::Coordinate lat_point(unit.lattice());
  const auto &basis = unit.basis();

  do {
    lat_point.frac() = grid_count().cast<double>();

    for (auto it = basis.begin(); it != basis.end(); ++it) {
      if (!site_filter(*it)) {
        continue;
      }

      xtal::Coordinate test(*it + lat_point);
      auto within_radius = [&](const xtal::Coordinate &coord) {
        return test.dist(coord) < max_radius;
      };
      if (std::any_of(basis.begin(), basis.end(), within_radius)) {
        *result++ = xtal::UnitCellCoord::from_coordinate(unit, test, xtal_tol);
      }
    }
  } while (++grid_count);
  return result;
}

/// \brief Output the neighborhood of sites within cutoff_radius of any sites in
/// the phenomenal
///
/// \param prim xtal::BasicStructure
/// \param phenomenal IntegralCluster
/// \param cutoff_radius The neighborhood distance cutoff
/// \param site_filter A filter function that returns true for UnitCellCoord
/// that
///        should be considered for the neighborhood
/// \param result Output iterator for container of UnitCellCoord
/// \param xtal_tol Crystallography tolerance used to contstruct UnitCellCoord
///
/// \returns Output iterator after generating the neighborhood
///
/// \ingroup IntegralCluster
///
template <typename OutputIterator>
OutputIterator neighborhood(xtal::BasicStructure const &prim,
                            IntegralCluster const &phenomenal,
                            double cutoff_radius,
                            SiteFilterFunction site_filter,
                            bool include_phenomenal_sites,
                            OutputIterator result, double xtal_tol) {
  int max_low_shift = 0;
  int max_high_shift = 0;
  for (auto it = phenomenal.begin(); it != phenomenal.end(); it++) {
    Eigen::Vector3l vec = it->unitcell();
    if (vec.maxCoeff() > max_high_shift) {
      max_high_shift = vec.maxCoeff();
    }
    if (vec.minCoeff() < max_low_shift) {
      max_low_shift = vec.minCoeff();
    }
  }

  // make grid counter
  auto dim = prim.lattice().enclose_sphere(cutoff_radius);
  Eigen::Vector3i ones(1, 1, 1);
  EigenCounter<Eigen::Vector3i> grid_count(
      -dim + (max_low_shift * ones).cast<int>(),
      dim + (max_high_shift * ones).cast<int>(), Eigen::Vector3i::Constant(1));

  /// lattice scaling
  xtal::Coordinate lat_point(prim.lattice());

  const auto &basis = prim.basis();

  do {
    lat_point.frac() = grid_count().cast<double>();

    for (auto it = basis.begin(); it != basis.end(); ++it) {
      if (!site_filter(*it)) {
        continue;
      }

      xtal::Coordinate test{*it + lat_point};

      if (!include_phenomenal_sites) {
        auto is_almost_equal = [&](const xtal::UnitCellCoord &uccoord) {
          return test.dist(uccoord.coordinate(prim)) < xtal_tol;
        };
        if (std::any_of(phenomenal.begin(), phenomenal.end(),
                        is_almost_equal)) {
          continue;
        }
      }

      auto within_radius = [&](const xtal::UnitCellCoord &uccoord) {
        return test.dist(uccoord.coordinate(prim)) < cutoff_radius;
      };
      if (std::any_of(phenomenal.begin(), phenomenal.end(), within_radius)) {
        *result++ = xtal::UnitCellCoord::from_coordinate(prim, test, xtal_tol);
      }
    }
  } while (++grid_count);
  return result;
}

}  // namespace

/// \brief Default constructor
///
/// Notes:
/// - Default site_filter is `dof_sites_filter` (include all sites with DoF)
ClusterSpecs::ClusterSpecs()
    : site_filter_method("dof_sites"), site_filter(dof_sites_filter()) {}

/// \brief Constructor
///
/// Notes:
/// - Default site_filter is `dof_sites_filter` (include all sites with DoF)
ClusterSpecs::ClusterSpecs(
    std::shared_ptr<xtal::BasicStructure const> const &_prim,
    std::shared_ptr<SymGroup const> const &_generating_group)
    : prim(_prim),
      generating_group(_generating_group),
      site_filter_method("dof_sites"),
      site_filter(dof_sites_filter()) {}

namespace ClusterSpecs_impl {

class DoFSitesFilter {
 public:
  DoFSitesFilter(std::vector<DoFKey> const &_dofs) : dofs(_dofs) {}

  bool operator()(xtal::Site const &site) {
    if (dofs.empty() &&
        (site.dof_size() != 0 || site.occupant_dof().size() > 1)) {
      return true;
    }
    for (DoFKey const &dof : dofs) {
      if (site.has_dof(dof)) {
        return true;
      } else if (dof == "occ" && site.occupant_dof().size() > 1) {
        return true;
      }
    }
    return false;
  }

  std::vector<DoFKey> dofs;
};

class AllClusters {
 public:
  bool operator()(ClusterInvariants const &invariants,
                  IntegralCluster const &clust) {
    return true;
  }
};

class MaxLengthClusterFilter {
 public:
  MaxLengthClusterFilter(double _max_length) : max_length(_max_length) {}

  bool operator()(ClusterInvariants const &invariants,
                  IntegralCluster const &clust) {
    if (clust.size() <= 1) {
      return true;
    }
    return invariants.distances().back() < max_length;
  }

 private:
  double max_length;
};

class EmptyNeighborhood {
 public:
  std::vector<xtal::UnitCellCoord> operator()(xtal::BasicStructure const &prim,
                                              SiteFilterFunction site_filter) {
    return std::vector<xtal::UnitCellCoord>{};
  }
};

class OriginNeighborhood {
 public:
  std::vector<xtal::UnitCellCoord> operator()(xtal::BasicStructure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    for (int i = 0; i < prim.basis().size(); ++i) {
      if (site_filter(prim.basis()[i])) {
        result.emplace_back(i, 0, 0, 0);
      }
    }
    return result;
  }
};

class MaxLengthNeighborhood {
 public:
  MaxLengthNeighborhood(double _max_length) : max_length(_max_length){};

  std::vector<xtal::UnitCellCoord> operator()(xtal::BasicStructure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    double xtal_tol = prim.lattice().tol();
    neighborhood(prim, max_length, site_filter, std::back_inserter(result),
                 xtal_tol);
    return result;
  }

 private:
  double max_length;
};

/// Generate a vector of UnitCellCoord that are within cutoff_radius distance to
/// any site in the phenomenal cluster
class CutoffRadiusNeighborhood {
 public:
  CutoffRadiusNeighborhood(IntegralCluster const &_phenomenal,
                           double _cutoff_radius,
                           bool _include_phenomenal_sites)
      : phenomenal(_phenomenal),
        cutoff_radius(_cutoff_radius),
        include_phenomenal_sites(_include_phenomenal_sites){};

  std::vector<xtal::UnitCellCoord> operator()(xtal::BasicStructure const &prim,
                                              SiteFilterFunction site_filter) {
    std::vector<xtal::UnitCellCoord> result;
    double xtal_tol = prim.lattice().tol();
    neighborhood(prim, phenomenal, cutoff_radius, site_filter,
                 include_phenomenal_sites, std::back_inserter(result),
                 xtal_tol);
    return result;
  }

 private:
  IntegralCluster phenomenal;
  double cutoff_radius;
  bool include_phenomenal_sites;
};

}  // namespace ClusterSpecs_impl

/// \brief Generate clusters using all Site
bool all_sites_filter(xtal::Site const &site) { return true; }

/// \brief Generate clusters using Site with site_occupant.size() > 1
bool alloy_sites_filter(xtal::Site const &site) {
  return site.occupant_dof().size() > 1;
}

/// \brief Generate clusters using Site with specified DoF
///
/// If dofs is empty, return true if Site has any continuous DoF or >1 allowed
/// occupant DoF If dofs is not empty, return true if Site has any of the DoF
/// types included. Use "occ" for / Site with >1 occupant allowed
SiteFilterFunction dof_sites_filter(std::vector<DoFKey> const &dofs) {
  return ClusterSpecs_impl::DoFSitesFilter{dofs};
}

/// Accept all clusters
ClusterFilterFunction all_clusters_filter() {
  return ClusterSpecs_impl::AllClusters{};
}

/// Accept clusters with max pair distance less than max_length
ClusterFilterFunction max_length_cluster_filter(double max_length) {
  return ClusterSpecs_impl::MaxLengthClusterFilter{max_length};
}

/// No sites (for null orbit, or global dof only)
CandidateSitesFunction empty_neighborhood() {
  return ClusterSpecs_impl::EmptyNeighborhood{};
}

/// Only sites in the origin unit cell {b, 0, 0, 0}
CandidateSitesFunction origin_neighborhood() {
  return ClusterSpecs_impl::OriginNeighborhood{};
}

/// Sites within max_length distance to any site in the origin unit cell {b, 0,
/// 0, 0}
CandidateSitesFunction max_length_neighborhood(double max_length) {
  return ClusterSpecs_impl::MaxLengthNeighborhood{max_length};
}

/// Sites within cutoff_radius distance to any site in the phenomenal cluster
CandidateSitesFunction cutoff_radius_neighborhood(
    IntegralCluster const &phenomenal, double cutoff_radius,
    bool include_phenomenal_sites) {
  return ClusterSpecs_impl::CutoffRadiusNeighborhood{phenomenal, cutoff_radius,
                                                     include_phenomenal_sites};
}

}  // namespace clust
}  // namespace CASM
