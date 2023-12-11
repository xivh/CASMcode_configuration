#include "casm/configuration/config_space_analysis.hh"

#include "casm/casm_io/Log.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "testdir.hh"
#include "teststructures.hh"

using namespace CASM;

class ConfigSpaceAnalysisTest : public testing::Test {
 protected:
  ConfigSpaceAnalysisTest() { log = Log(std::cout, Log::standard, true); }

  void make_prim(xtal::BasicStructure const &_xtal_prim) {
    xtal_prim = std::make_shared<xtal::BasicStructure>(_xtal_prim);
    prim = std::make_shared<config::Prim>(xtal_prim);
  }

  void read_prim_file(std::string file_name) {
    xtal_prim = std::make_shared<xtal::BasicStructure>(
        read_prim(test::data_file("configuration", file_name), TOL));
    prim = std::make_shared<config::Prim>(xtal_prim);
  }

  Index to_site_index(std::shared_ptr<config::Supercell const> const &supercell,
                      xtal::UnitCellCoord unitcellcoord) {
    return supercell->unitcellcoord_index_converter(unitcellcoord);
  }

  int &occ(config::Configuration &config, Index site_index) {
    return config.dof_values.occupation(site_index);
  }

  int &occ(config::Configuration &config, xtal::UnitCellCoord unitcellcoord) {
    return occ(config, to_site_index(config.supercell, unitcellcoord));
  }

  void build_configurations_1() {
    {
      Eigen::Matrix3l T;
      T << 1, 0, 0,  //
          0, 1, 0,   //
          0, 0, 1;   //
      auto shared_supercell =
          std::make_shared<config::Supercell const>(prim, T);
      config::Configuration config(shared_supercell);

      occ(config, 0) = 0;
      configurations.emplace("A1", config);

      occ(config, 0) = 1;
      configurations.emplace("B1", config);
    }

    {
      Eigen::Matrix3l T;
      T << -1, 1, 1,  //
          1, -1, 1,   //
          1, 1, -1;   //
      auto shared_supercell =
          std::make_shared<config::Supercell const>(prim, T);

      config::Configuration config(shared_supercell);
      Eigen::VectorXi occ(4);

      occ << 1, 0, 0, 0;
      config.dof_values.occupation = occ;
      configurations.emplace("A3B1", config);

      occ << 0, 1, 1, 1;
      config.dof_values.occupation = occ;
      configurations.emplace("A1B3", config);
    }
  }

  void expect_same_basis_vectors(Eigen::MatrixXd const &basis,
                                 Eigen::MatrixXd const &expected) {
    EXPECT_EQ(basis.cols(), expected.cols());
    bool all_basis_vectors_found = true;
    for (Index c_basis = 0; c_basis < basis.cols(); ++c_basis) {
      bool found = false;
      for (Index c_expected = 0; c_expected < expected.cols(); ++c_expected) {
        if (almost_equal(basis.col(c_basis), expected.col(c_expected))) {
          found = true;
          break;
        }
      }
      if (!found) {
        all_basis_vectors_found = false;
        break;
      }
    }
    EXPECT_TRUE(all_basis_vectors_found) << "basis: \n"
                                         << basis << "\nexpected: \n"
                                         << expected << std::endl;
  }

  std::shared_ptr<xtal::BasicStructure const> xtal_prim;
  std::shared_ptr<config::Prim const> prim;

  // Config space analysis options:

  std::map<std::string, config::Configuration> configurations;
  std::optional<std::vector<DoFKey>> dofs = std::nullopt;
  std::optional<bool> exclude_homogeneous_modes = std::nullopt;
  bool include_default_occ_modes = false;
  std::optional<std::map<int, int>> sublattice_index_to_default_occ =
      std::nullopt;
  std::optional<std::map<Index, int>> site_index_to_default_occ = std::nullopt;
  double tol = TOL;

  /// Logging
  std::optional<Log> log;
};

TEST_F(ConfigSpaceAnalysisTest, Test1) {
  make_prim(test::FCC_binary_prim());
  build_configurations_1();

  // Perform config space analysis
  std::map<DoFKey, config::ConfigSpaceAnalysisResults> results =
      config::config_space_analysis(
          configurations, dofs, exclude_homogeneous_modes,
          include_default_occ_modes, sublattice_index_to_default_occ,
          site_index_to_default_occ, tol);

  // Check equivalents
  auto const &eq = results.at("occ").equivalent_configurations;
  EXPECT_EQ(eq.at("A1").size(), 1);
  EXPECT_EQ(eq.at("B1").size(), 1);
  EXPECT_EQ(eq.at("A1B3").size(), 4);
  EXPECT_EQ(eq.at("A3B1").size(), 4);

  // Check basis
  Eigen::MatrixXd const &basis =
      results.at("occ").symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);

  Eigen::MatrixXd expected(8, 4);
  expected.col(0) << 0.0, 0.0, 0.0, sqrt(6.) / 3., 0.0, -sqrt(6.) / 6., 0.0,
      -sqrt(6.) / 6.;
  expected.col(1) << 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0,
      sqrt(2.0) / 2.0;
  expected.col(2) << 0.0, -sqrt(3.) / 2., 0.0, sqrt(3.) / 6., 0.0,
      sqrt(3.) / 6., 0.0, sqrt(3.) / 6.;
  expected.col(3) << 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5;
  expect_same_basis_vectors(basis, expected);
}

TEST_F(ConfigSpaceAnalysisTest, Test2) {
  make_prim(test::FCC_binary_prim());
  build_configurations_1();
  sublattice_index_to_default_occ = std::map<int, int>({{0, 1}});

  // Perform config space analysis
  std::map<DoFKey, config::ConfigSpaceAnalysisResults> results =
      config::config_space_analysis(
          configurations, dofs, exclude_homogeneous_modes,
          include_default_occ_modes, sublattice_index_to_default_occ,
          site_index_to_default_occ, tol);

  // Check equivalents
  auto const &eq = results.at("occ").equivalent_configurations;
  EXPECT_EQ(eq.at("A1").size(), 1);
  EXPECT_EQ(eq.at("B1").size(), 1);
  EXPECT_EQ(eq.at("A1B3").size(), 4);
  EXPECT_EQ(eq.at("A3B1").size(), 4);

  // Check basis
  Eigen::MatrixXd const &basis =
      results.at("occ").symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);

  Eigen::MatrixXd expected(8, 4);
  expected.col(0) << 0.0, 0.0, sqrt(6.) / 3., 0.0, -sqrt(6.) / 6., 0.0,
      -sqrt(6.) / 6., 0.0;
  expected.col(1) << 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0,
      0.0;
  expected.col(2) << -sqrt(3.) / 2., 0.0, sqrt(3.) / 6., 0.0, sqrt(3.) / 6.,
      0.0, sqrt(3.) / 6., 0.0;
  expected.col(3) << 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0;
  expect_same_basis_vectors(basis, expected);
}

TEST_F(ConfigSpaceAnalysisTest, Test3) {
  make_prim(test::FCC_binary_prim());
  build_configurations_1();
  site_index_to_default_occ =
      std::map<Index, int>({{0, 1}, {1, 1}, {2, 0}, {3, 0}});

  // Perform config space analysis
  std::map<DoFKey, config::ConfigSpaceAnalysisResults> results =
      config::config_space_analysis(
          configurations, dofs, exclude_homogeneous_modes,
          include_default_occ_modes, sublattice_index_to_default_occ,
          site_index_to_default_occ, tol);

  // Check equivalents
  auto const &eq = results.at("occ").equivalent_configurations;
  EXPECT_EQ(eq.at("A1").size(), 1);
  EXPECT_EQ(eq.at("B1").size(), 1);
  EXPECT_EQ(eq.at("A1B3").size(), 4);
  EXPECT_EQ(eq.at("A3B1").size(), 4);

  // Check basis
  Eigen::MatrixXd const &basis =
      results.at("occ").symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);

  Eigen::MatrixXd expected(8, 4);
  expected.col(0) << 0.0, 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0,
      sqrt(2.0) / 2.0;
  expected.col(1) << -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0, 0.0,
      0.0;
  expected.col(2) << -0.5, 0.0, -0.5, 0.0, 0.0, 0.5, 0.0, 0.5;
  expected.col(3) << 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5;
  EXPECT_TRUE(almost_equal(basis, expected));
}

TEST_F(ConfigSpaceAnalysisTest, Test4) {
  make_prim(test::FCC_binary_prim());
  build_configurations_1();
  site_index_to_default_occ =
      std::map<Index, int>({{0, 0}, {1, 0}, {2, 1}, {3, 1}});

  // Perform config space analysis
  std::map<DoFKey, config::ConfigSpaceAnalysisResults> results =
      config::config_space_analysis(
          configurations, dofs, exclude_homogeneous_modes,
          include_default_occ_modes, sublattice_index_to_default_occ,
          site_index_to_default_occ, tol);

  // Check equivalents
  auto const &eq = results.at("occ").equivalent_configurations;
  EXPECT_EQ(eq.at("A1").size(), 1);
  EXPECT_EQ(eq.at("B1").size(), 1);
  EXPECT_EQ(eq.at("A1B3").size(), 4);
  EXPECT_EQ(eq.at("A3B1").size(), 4);

  // Check basis
  Eigen::MatrixXd const &basis =
      results.at("occ").symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);

  Eigen::MatrixXd expected(8, 4);
  expected.col(0) << 0.0, 0.0, 0.0, 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0,
      0.0;
  expected.col(1) << 0.0, -sqrt(2.0) / 2.0, 0.0, sqrt(2.0) / 2.0, 0.0, 0.0, 0.0,
      0.0;
  expected.col(2) << 0.0, -0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.0;
  expected.col(3) << 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0;
  EXPECT_TRUE(almost_equal(basis, expected));
}