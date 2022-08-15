#include "casm/configuration/clusterography/orbits.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/ClusterInvariants.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
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

void print(std::shared_ptr<clust::SymGroup const> const &factor_group,
           xtal::Lattice const &lattice) {
  xtal::SymInfoOptions opt{FRAC};
  std::cout << "factor_group: " << std::endl;
  for (auto const &op : factor_group->element) {
    xtal::SymInfo syminfo{op, lattice};
    std::cout << to_brief_unicode(syminfo, opt) << std::endl;
  }
  std::cout << std::endl;
}

void print(std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
           xtal::BasicStructure const &prim) {
  jsonParser json;
  json.put_array();
  for (auto const &orbit : orbits) {
    jsonParser tjson;
    tjson.put_array();
    for (auto const &cluster : orbit) {
      jsonParser cjson;
      to_json(cluster, cjson, prim);
      tjson.push_back(cjson);
    }
    json.push_back(tjson);
  }
  std::cout << json << std::endl;
}

void print_prototypes(
    std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
    xtal::BasicStructure const &prim) {
  jsonParser json;
  json.put_array();
  for (auto const &orbit : orbits) {
    jsonParser tjson;
    to_json(orbit.size(), tjson["orbit_size"]);
    to_json(*orbit.begin(), tjson["prototype"], prim);
    json.push_back(tjson);
  }
  std::cout << json << std::endl;
  std::cout << "orbits.size(): " << orbits.size() << std::endl;
  std::cout << "cluster_size: ";
  for (auto const &orbit : orbits) {
    std::cout << orbit.begin()->size() << ", ";
  }
  std::cout << std::endl;
  std::cout << "orbit_size: ";
  for (auto const &orbit : orbits) {
    std::cout << orbit.size() << ", ";
  }
  std::cout << std::endl;
}

void print_local_prototypes(
    std::vector<std::set<CASM::clust::IntegralCluster>> const &orbits,
    xtal::BasicStructure const &prim,
    clust::IntegralCluster const &phenomenal) {
  jsonParser json;
  json.put_array();
  for (auto const &orbit : orbits) {
    jsonParser tjson;
    to_json(orbit.size(), tjson["orbit_size"]);
    to_json(*orbit.begin(), tjson["prototype"], prim);

    clust::ClusterInvariants invariants(*orbit.begin(), phenomenal, prim);
    if (invariants.phenomenal_distances().size()) {
      tjson["phenomenal_to_cluster"]["min_length"] =
          invariants.phenomenal_distances().front();
      tjson["phenomenal_to_cluster"]["max_length"] =
          invariants.phenomenal_distances().back();
    } else {
      tjson["phenomenal_to_cluster"]["min_length"] = 0.0;
      tjson["phenomenal_to_cluster"]["max_length"] = 0.0;
    }
    tjson["phenomenal_to_cluster"]["distances"] =
        invariants.phenomenal_distances();

    json.push_back(tjson);
  }
  std::cout << json << std::endl;
  std::cout << "orbits.size(): " << orbits.size() << std::endl;
  std::cout << "cluster_size: ";
  for (auto const &orbit : orbits) {
    std::cout << orbit.begin()->size() << ", ";
  }
  std::cout << std::endl;
  std::cout << "orbit_size: ";
  for (auto const &orbit : orbits) {
    std::cout << orbit.size() << ", ";
  }
  std::cout << std::endl;
}

void print_phenomenal(clust::IntegralCluster const &cluster) {
  std::cout << "phenomenal:" << std::endl;
  for (auto const &site : cluster.elements()) {
    std::cout << "{" << site.sublattice() << "," << site.unitcell()(0) << ","
              << site.unitcell()(1) << "," << site.unitcell()(2) << "}"
              << std::endl;
  }
}

void print_neighborhood(std::set<xtal::UnitCellCoord> const &neighborhood) {
  std::cout << "neighborhood:" << std::endl;
  for (auto const &site : neighborhood) {
    std::cout << "{" << site.sublattice() << "," << site.unitcell()(0) << ","
              << site.unitcell()(1) << "," << site.unitcell()(2) << "}"
              << std::endl;
  }
}

}  // namespace test

// test FCC_binary_prim
TEST(PrimPeriodicOrbitTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  auto factor_group = sym_info::make_factor_group(*prim);
  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
  std::vector<double> max_length = {0, 0, 4.01, 4.01};
  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

  auto orbits =
      make_prim_periodic_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                max_length, custom_generators);
  // test::print(orbits, *prim);
  // test::print_prototypes(orbits, *prim);
  EXPECT_EQ(orbits.size(), 6);

  std::vector<Index> cluster_size = {0, 1, 2, 2, 3, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }
}

// test ZrO (some sites with no DoF)
TEST(PrimPeriodicOrbitTest, Test2) {
  auto prim = std::make_shared<xtal::BasicStructure const>(test::ZrO_prim());
  auto factor_group = sym_info::make_factor_group(*prim);
  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
  std::vector<double> max_length = {0, 0, 5.17, 5.17};
  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

  auto orbits =
      make_prim_periodic_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                max_length, custom_generators);
  // test::print_prototypes(orbits, *prim);
  EXPECT_EQ(orbits.size(), 12);

  std::vector<Index> cluster_size = {0, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }
}

// test ZrO w/all_sites_filter
TEST(PrimPeriodicOrbitTest, Test3) {
  auto prim = std::make_shared<xtal::BasicStructure const>(test::ZrO_prim());
  auto factor_group = sym_info::make_factor_group(*prim);
  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::all_sites_filter;
  std::vector<double> max_length = {0, 0, 5.17, 5.17};
  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

  auto orbits =
      make_prim_periodic_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                max_length, custom_generators);
  // test::print_prototypes(orbits, *prim);
  EXPECT_EQ(orbits.size(), 64);

  std::vector<Index> cluster_size = {
      0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }
}

// test FCC_binary prim w/custom clusters
TEST(PrimPeriodicOrbitTest, Test4) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  auto factor_group = sym_info::make_factor_group(*prim);
  auto unitcellcoord_symgroup_rep =
      sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
  clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
  std::vector<double> max_length = {};
  std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {
      clust::IntegralClusterOrbitGenerator(clust::IntegralCluster({
                                               xtal::UnitCellCoord(0, 0, 0, 0),
                                               xtal::UnitCellCoord(0, 0, 0, 1),
                                               xtal::UnitCellCoord(0, 0, 1, 0),
                                           }),
                                           true  // include_subclusters
                                           )};

  auto orbits =
      make_prim_periodic_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                max_length, custom_generators);
  // test::print(orbits, *prim);
  // test::print_prototypes(orbits, *prim);
  EXPECT_EQ(orbits.size(), 4);

  std::vector<Index> cluster_size = {0, 1, 2, 3};
  auto cluster_size_it = cluster_size.begin();
  for (auto const &orbit : orbits) {
    EXPECT_EQ(orbit.begin()->size(), *cluster_size_it++);
  }
}
