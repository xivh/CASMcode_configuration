#include "casm/configuration/enumeration/perturbations.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

using namespace CASM;

class FCCBinaryPerturbationsTest : public testing::Test {
 protected:
  FCCBinaryPerturbationsTest() {
    auto basicstructure =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    prim = std::make_shared<config::Prim>(basicstructure);
  }

  /// \brief Make configurations that are distinct occupation perturbations of
  /// the given clusters
  std::set<config::Configuration> make_distinct_perturbations(
      config::Configuration const &configuration,
      std::vector<clust::IntegralCluster> const &clusters) {
    auto const &supercell = configuration.supercell;
    auto const &prim = supercell->prim;
    auto const &xtal_prim = prim->basicstructure;
    std::vector<std::set<clust::IntegralCluster>> orbits;
    for (auto const &cluster : clusters) {
      orbits.emplace_back(make_prim_periodic_orbit(
          cluster, prim->sym_info.unitcellcoord_symgroup_rep));
    }

    auto orbits_as_indices = clust::make_orbits_as_indices(
        orbits, supercell->unitcellcoord_index_converter);
    {
      jsonParser json;
      to_json(orbits_as_indices, json["orbits_as_indices"]);
      std::cout << json << std::endl;
    }

    auto distinct_cluster_sites =
        make_distinct_cluster_sites(configuration, orbits_as_indices);
    {
      jsonParser json;
      to_json(distinct_cluster_sites, json["distinct_cluster_sites"]);
      std::cout << json << std::endl;
    }
    return config::make_distinct_perturbations(configuration,
                                               distinct_cluster_sites);
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(FCCBinaryPerturbationsTest, Test1) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}, {0, 1, 0, 0}})});
  config::Configuration configuration(supercell);
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 3);
}

TEST_F(FCCBinaryPerturbationsTest, Test2) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * m^3 total sites
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  L *= 3;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}, {0, 1, 0, 0}})});
  config::Configuration configuration(supercell);
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 3);
}

TEST_F(FCCBinaryPerturbationsTest, Test3) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * 3
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 12.;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}, {0, 1, 0, 0}})});
  config::Configuration configuration(supercell);
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 4);
}

TEST_F(FCCBinaryPerturbationsTest, Test4) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * 6
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 8., 0.;
  L.col(2) << 0., 0., 12.;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}, {0, 1, 0, 0}})});
  config::Configuration configuration(supercell);
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 5);
}

TEST_F(FCCBinaryPerturbationsTest, Test5) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * m^3 total sites
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  L *= 1;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}})});
  config::Configuration configuration(supercell);
  configuration.dof_values.occupation(0) = 1;
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 3);
}

TEST_F(FCCBinaryPerturbationsTest, Test6) {
  using namespace clust;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * 2 total sites
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 8.;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  std::vector<clust::IntegralCluster> clusters(
      {clust::IntegralCluster({{0, 0, 0, 0}})});
  config::Configuration configuration(supercell);
  configuration.dof_values.occupation(0) = 1;
  std::set<config::Configuration> perturbations =
      make_distinct_perturbations(configuration, clusters);
  for (auto const &c : perturbations) {
    std::cout << c.dof_values.occupation.transpose() << std::endl;
  }
  EXPECT_EQ(perturbations.size(), 6);
}

TEST_F(FCCBinaryPerturbationsTest, Test7) {
  using namespace clust;

  Eigen::Matrix3d motif_L;
  // conventional 4-atom fcc supercell * 2 total sites
  motif_L.col(0) << 4., 0., 0.;
  motif_L.col(1) << 0., 4., 0.;
  motif_L.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<config::Supercell const>(prim, xtal::Lattice(motif_L));

  // L12 config
  config::Configuration motif(motif_supercell);
  motif.dof_values.occupation(0) = 1;

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell * 2 total sites
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  L *= 2;
  supercell = std::make_shared<config::Supercell const>(prim, xtal::Lattice(L));

  // L12 super-config
  config::Configuration configuration =
      config::copy_configuration(motif, supercell);

  // point-cluster perturbations
  {
    std::vector<clust::IntegralCluster> clusters(
        {clust::IntegralCluster({{0, 0, 0, 0}})});

    std::set<config::Configuration> perturbations =
        make_distinct_perturbations(configuration, clusters);
    for (auto const &c : perturbations) {
      std::cout << c.dof_values.occupation.transpose() << std::endl;
    }
    EXPECT_EQ(perturbations.size(), 3);
  }

  // 1NN pair-cluster perturbations
  {
    std::vector<clust::IntegralCluster> clusters(
        {clust::IntegralCluster({{0, 0, 0, 0}, {0, 1, 0, 0}})});

    std::set<config::Configuration> perturbations =
        make_distinct_perturbations(configuration, clusters);
    for (auto const &c : perturbations) {
      std::cout << c.dof_values.occupation.transpose() << std::endl;
    }
    EXPECT_EQ(perturbations.size(), 5);
  }
}