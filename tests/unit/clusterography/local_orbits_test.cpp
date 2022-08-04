#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug
#include "casm/casm_io/container/stream_io.hh"
#include "casm/crystallography/SymInfo.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"

using namespace CASM;

namespace test {

void print(std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
           xtal::BasicStructure const &prim);

void print_prototypes(
    std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
    xtal::BasicStructure const &prim);

void print_local_prototypes(
    std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
    xtal::BasicStructure const &prim, clust::IntegralCluster const &phenomenal);

}  // namespace test

// test FCC_binary_prim - phenomenal == point cluster
TEST(LocalOrbitTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  auto factor_group = sym_info::make_factor_group(*prim);
  auto factor_group_unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::IntegralCluster phenomenal({xtal::UnitCellCoord(0, 0, 0, 0)});
  auto cluster_group = make_cluster_group(
      phenomenal, factor_group, prim->lattice().lat_column_mat(),
      factor_group_unitcellcoord_symgroup_rep);
  EXPECT_EQ(cluster_group->element.size(), 48);

  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(cluster_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
  std::vector<double> max_length = {0, 0, 4.01, 4.01};
  std::vector<double> cutoff_radius = {0, 4.01, 4.01, 4.01};
  bool include_phenomenal_sites = false;

  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

  auto orbits = make_local_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                  max_length, custom_generators, phenomenal,
                                  cutoff_radius, include_phenomenal_sites);
  // test::print(orbits, *prim);
  // test::print_local_prototypes(orbits, *prim, phenomenal);
  EXPECT_EQ(orbits.size(), 10);

  std::vector<Index> cluster_size = {0, 1, 1, 2, 2, 2, 3, 3, 3, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }

  std::vector<Index> orbit_size = {1, 12, 6, 24, 24, 12, 8, 24, 24, 12};
  auto orbit_size_it = orbit_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.size(), *orbit_size_it++);
  }
}

// test FCC_binary_prim - phenomenal == 1NN pair cluster
TEST(LocalOrbitTest, Test2) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  auto factor_group = sym_info::make_factor_group(*prim);

  auto factor_group_unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::IntegralCluster phenomenal(
      {xtal::UnitCellCoord(0, 0, 0, 0), xtal::UnitCellCoord(0, 0, 1, 0)});
  auto cluster_group = make_cluster_group(
      phenomenal, factor_group, prim->lattice().lat_column_mat(),
      factor_group_unitcellcoord_symgroup_rep);
  EXPECT_EQ(cluster_group->element.size(), 8);

  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(cluster_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
  std::vector<double> max_length = {0, 0, 4.01, 4.01};
  std::vector<double> cutoff_radius = {0, 4.01, 4.01, 4.01};
  bool include_phenomenal_sites = false;

  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

  auto orbits = make_local_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                  max_length, custom_generators, phenomenal,
                                  cutoff_radius, include_phenomenal_sites);
  // test::print(orbits, *prim);
  // test::print_local_prototypes(orbits, *prim, phenomenal);
  EXPECT_EQ(orbits.size(), 42);

  std::vector<Index> cluster_size = {0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
                                     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
                                     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }

  std::vector<Index> orbit_size = {1, 4, 4, 8, 4, 2, 4, 2, 8, 2, 8, 8, 8, 4,
                                   8, 2, 8, 4, 8, 4, 2, 8, 4, 4, 4, 4, 8, 8,
                                   4, 4, 4, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 4};
  auto orbit_size_it = orbit_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.size(), *orbit_size_it++);
  }
}
