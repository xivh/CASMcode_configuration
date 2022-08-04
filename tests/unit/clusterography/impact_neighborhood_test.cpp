#include "casm/configuration/clusterography/impact_neighborhood.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

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

void print_neighborhood(std::set<xtal::UnitCellCoord> const &neighborhood);

void print_phenomenal(clust::IntegralCluster const &cluster);

}  // namespace test

// test FCC_binary_prim
class FlowerNeighborhoodTest : public testing::Test {
 protected:
  std::vector<std::set<clust::IntegralCluster>> orbits;

  FlowerNeighborhoodTest() {
    auto prim =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    auto factor_group = sym_info::make_factor_group(*prim);
    auto unitcellcoord_symgroup_rep =
        sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);
    clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
    std::vector<double> max_length = {0, 0, 4.01, 4.01};
    std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

    orbits =
        make_prim_periodic_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                                  max_length, custom_generators);
  }
};

TEST_F(FlowerNeighborhoodTest, Test1) {
  // 1-site phenomenal
  clust::IntegralCluster phenomenal({xtal::UnitCellCoord(0, 0, 0, 0)});

  // neighborhood
  std::set<xtal::UnitCellCoord> neighborhood;
  add_to_flower_neighborhood(phenomenal, neighborhood, orbits);

  // test::print_neighborhood(neighborhood);

  // 1 point + 12 1NN + 6 2NN = 19 (phenomenal site included via point cluster)
  EXPECT_EQ(neighborhood.size(), 19);
}

TEST_F(FlowerNeighborhoodTest, Test2) {
  // 2-site phenomenal (1NN)
  clust::IntegralCluster phenomenal(
      {xtal::UnitCellCoord(0, 0, 0, 0), xtal::UnitCellCoord(0, 1, 0, 0)});

  // neighborhood
  std::set<xtal::UnitCellCoord> neighborhood;
  add_to_flower_neighborhood(phenomenal, neighborhood, orbits);

  // test::print_neighborhood(neighborhood);

  // 19 + 19 - 10 overlap
  EXPECT_EQ(neighborhood.size(), 28);
}

// test FCC_binary_prim - phenomenal == point cluster
class LocalNeighborhoodTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<clust::SymGroup const> factor_group;
  std::vector<xtal::UnitCellCoordRep> factor_group_unitcellcoord_symgroup_rep;
  clust::IntegralCluster phenomenal;
  std::vector<std::set<clust::IntegralCluster>> orbits;

  LocalNeighborhoodTest() {
    // prim & factor group
    prim =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    factor_group = sym_info::make_factor_group(*prim);
    factor_group_unitcellcoord_symgroup_rep =
        sym_info::make_unitcellcoord_symgroup_rep(factor_group->element, *prim);

    // phenomenal cluster & cluster group
    phenomenal = clust::IntegralCluster(
        {xtal::UnitCellCoord(0, 0, 0, 0), xtal::UnitCellCoord(0, 1, 0, 0)});
    auto cluster_group = make_cluster_group(
        phenomenal, factor_group, prim->lattice().lat_column_mat(),
        factor_group_unitcellcoord_symgroup_rep);
    EXPECT_EQ(cluster_group->element.size(), 8);
    auto unitcellcoord_symgroup_rep = sym_info::make_unitcellcoord_symgroup_rep(
        cluster_group->element, *prim);

    // local orbit parameters
    clust::SiteFilterFunction site_filter = clust::dof_sites_filter();
    std::vector<double> max_length = {0, 0, 4.01, 4.01};
    std::vector<double> cutoff_radius = {0, 4.01, 4.01, 4.01};
    bool include_phenomenal_sites = false;
    std::vector<clust::IntegralClusterOrbitGenerator> custom_generators = {};

    orbits = make_local_orbits(prim, unitcellcoord_symgroup_rep, site_filter,
                               max_length, custom_generators, phenomenal,
                               cutoff_radius, include_phenomenal_sites);
  }
};

TEST_F(LocalNeighborhoodTest, Test1) {
  // neighborhood
  std::set<xtal::UnitCellCoord> neighborhood;
  add_to_local_neighborhood(neighborhood, orbits);

  // test::print_neighborhood(neighborhood);

  // 19 + 19 - 10 overlap - 2 (does not include phenomenal sites)
  EXPECT_EQ(neighborhood.size(), 26);
}

TEST_F(LocalNeighborhoodTest, Test2) {
  // neighborhood
  std::set<xtal::UnitCellCoord> neighborhood;
  add_to_local_neighborhood(neighborhood, orbits);

  // get equivalent phenomenal clusters & neighborhoods...

  // get phenomenal cluster orbit
  auto const &symgroup_rep = factor_group_unitcellcoord_symgroup_rep;
  auto phenom_orbit = make_prim_periodic_orbit(phenomenal, symgroup_rep);

  /// The indices equivalence_map[i] are
  ///     the indices of the group elements that transform the first
  ///     element in the orbit into the i-th element in the orbit.
  std::vector<std::vector<Index>> equivalence_map = group::make_equivalence_map(
      phenom_orbit, symgroup_rep.begin(), symgroup_rep.end(),
      clust::prim_periodic_integral_cluster_copy_apply);

  /// Generate all equivalent phenomenal & neighborhoods in the origin cell
  std::vector<clust::IntegralCluster> equiv_phenomenal;
  std::vector<std::set<xtal::UnitCellCoord>> equiv_neighborhood;
  for (Index i = 0; i < equivalence_map.size(); ++i) {
    xtal::UnitCellCoordRep const &rep = symgroup_rep[equivalence_map[i][0]];

    // transform phenomal:
    //
    //     phenomenal' = op*phenomenal + trans,
    //
    // where `trans` translates phenomenal' so the first site is
    // in the origin cell
    clust::IntegralCluster tphenomenal = copy_apply(rep, phenomenal);
    tphenomenal.sort();
    xtal::UnitCell trans = -tphenomenal[0].unitcell();
    tphenomenal += trans;

    // transform neighborhood sites:
    // site' = op*site + trans
    std::set<xtal::UnitCellCoord> tneighborhood;
    for (auto const &site : neighborhood) {
      tneighborhood.emplace(copy_apply(rep, site) + trans);
    }

    // std::cout << "--- i: " << i << " ---" << std::endl;
    // test::print_phenomenal(tphenomenal);
    // test::print_neighborhood(tneighborhood);

    equiv_phenomenal.emplace_back(std::move(tphenomenal));
    equiv_neighborhood.emplace_back(std::move(tneighborhood));
  }

  EXPECT_EQ(equiv_phenomenal.size(), 6);
  EXPECT_EQ(equiv_neighborhood.size(), 6);
}
