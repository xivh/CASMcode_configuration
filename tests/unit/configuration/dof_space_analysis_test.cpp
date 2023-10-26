#include "casm/configuration/dof_space_analysis.hh"

#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "testdir.hh"
#include "teststructures.hh"

using namespace CASM;

class DoFSpaceAnalysisTest : public testing::Test {
 protected:
  DoFSpaceAnalysisTest() { log = Log(std::cout, Log::standard, true); }

  void make_prim(xtal::BasicStructure const &_xtal_prim) {
    xtal_prim = std::make_shared<xtal::BasicStructure>(_xtal_prim);
    prim = std::make_shared<config::Prim>(xtal_prim);
  }

  void read_prim_file(std::string file_name) {
    xtal_prim = std::make_shared<xtal::BasicStructure>(
        read_prim(test::data_file("configuration", file_name), TOL));
    prim = std::make_shared<config::Prim>(xtal_prim);
  }

  void make_prim_dof_space(DoFKey dof_key) {
    Eigen::Matrix3l T;
    T << 1, 0, 0,  //
        0, 1, 0,   //
        0, 0, 1;   //
    transformation_matrix_to_super = T;
    dof_space = std::make_unique<clexulator::DoFSpace>(
        dof_key, xtal_prim, transformation_matrix_to_super, sites, basis);
  }

  void make_dof_space(DoFKey dof_key) {
    dof_space = std::make_unique<clexulator::DoFSpace>(
        dof_key, xtal_prim, transformation_matrix_to_super, sites, basis);
  }

  std::shared_ptr<xtal::BasicStructure const> xtal_prim;
  std::shared_ptr<config::Prim const> prim;

  // Initial DoFSpace construction:

  std::optional<Eigen::Matrix3l> transformation_matrix_to_super;
  std::optional<std::set<Index>> sites;
  std::optional<Eigen::MatrixXd> basis;
  std::unique_ptr<clexulator::DoFSpace> dof_space;

  // DoF Space analysis options:

  /// Specify symmetry via configuration factor group
  std::optional<config::Configuration> configuration;

  /// Exclude homogeneous modes. If this is null (default),
  /// exclude homogeneous modes for dof==\"disp\" only.
  std::optional<bool> exclude_homogeneous_modes;

  /// Include the dof component for the default occupation value on
  /// each site with occupation DoF. The default is to exclude these
  /// modes because they are not independent. This parameter is only
  /// checked dof==\"occ\".
  bool include_default_occ_modes = false;

  std::optional<std::map<int, int>> sublattice_index_to_default_occ =
      std::nullopt;

  std::optional<std::map<Index, int>> site_index_to_default_occ = std::nullopt;

  /// If true, calculate the irreducible wedges for the vector space.
  /// This may take a long time.
  bool calc_wedges = false;

  /// Logging
  std::optional<Log> log;
};

TEST_F(DoFSpaceAnalysisTest, Test1) {
  make_prim(test::FCC_binary_prim());
  make_prim_dof_space("occ");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  //  std::vector<Eigen::MatrixXd> const &symgroup_rep =
  //      symmetry_report.symgroup_rep;
  //  std::vector<irreps::SubWedge> const &irreducible_wedge =
  //      symmetry_report.irreducible_wedge;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;
  //  std::vector<std::string> const &axis_glossary =
  //  symmetry_report.axis_glossary;

  EXPECT_EQ(irreps.size(), 1);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 2);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 1);
}

TEST_F(DoFSpaceAnalysisTest, Test2) {
  make_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1,  //
      1, -1, 1,   //
      1, 1, -1;   //
  transformation_matrix_to_super = T;
  make_dof_space("occ");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 2);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 8);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 4);

  // Check basis
  Eigen::MatrixXd const &basis = results.symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);
}

TEST_F(DoFSpaceAnalysisTest, Test2a) {
  make_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1,  //
      1, -1, 1,   //
      1, 1, -1;   //
  transformation_matrix_to_super = T;
  make_dof_space("occ");

  site_index_to_default_occ =
      std::map<Index, int>({{0, 0}, {1, 0}, {2, 0}, {3, 1}});

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 4);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 8);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 8);

  // Check basis
  Eigen::MatrixXd const &basis = results.symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 8);
}

TEST_F(DoFSpaceAnalysisTest, Test2b) {
  make_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1,  //
      1, -1, 1,   //
      1, 1, -1;   //
  transformation_matrix_to_super = T;
  make_dof_space("occ");

  sublattice_index_to_default_occ = std::map<int, int>({{0, 1}});

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 2);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 8);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 4);

  // Check basis
  Eigen::MatrixXd const &basis = results.symmetry_adapted_dof_space.basis;
  //  std::cout << "basis:" << std::endl;
  //  std::cout << basis << std::endl;
  EXPECT_EQ(basis.rows(), 8);
  EXPECT_EQ(basis.cols(), 4);
}

TEST_F(DoFSpaceAnalysisTest, Test3) {
  read_prim_file("prim_test_occ_fix_corner.json");
  make_prim_dof_space("occ");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 2);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 7);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 3);
}

TEST_F(DoFSpaceAnalysisTest, Test4) {
  make_prim(test::FCC_binary_disp_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1,  //
      1, -1, 1,   //
      1, 1, -1;   //
  transformation_matrix_to_super = T;
  make_dof_space("disp");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 2);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 12);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 9);
}

TEST_F(DoFSpaceAnalysisTest, Test5) {
  read_prim_file("prim_ABC2.json");
  make_prim_dof_space("disp");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 15);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 18);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 15);
}

TEST_F(DoFSpaceAnalysisTest, Test6) {
  read_prim_file("prim_test_disp_fix_corner.json");
  make_prim_dof_space("disp");

  // Perform DoF space analysis
  config::DoFSpaceAnalysisResults results = config::dof_space_analysis(
      *dof_space, prim, configuration, exclude_homogeneous_modes,
      include_default_occ_modes, sublattice_index_to_default_occ,
      site_index_to_default_occ, calc_wedges, log);

  // Check results
  irreps::VectorSpaceSymReport const &symmetry_report = results.symmetry_report;
  std::vector<irreps::IrrepInfo> const &irreps = symmetry_report.irreps;
  Eigen::MatrixXd const &symmetry_adapted_subspace =
      symmetry_report.symmetry_adapted_subspace;

  EXPECT_EQ(irreps.size(), 3);
  EXPECT_EQ(symmetry_adapted_subspace.rows(), 9);
  EXPECT_EQ(symmetry_adapted_subspace.cols(), 9);
}