#include "casm/configuration/canonical_form.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class CanonicalFormFCCTest : public testing::Test {
 protected:
  CanonicalFormFCCTest() {
    prim = config::make_shared_prim(test::FCC_binary_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(CanonicalFormFCCTest, TestSupercell1) {
  Eigen::Matrix3d S;
  S << 4., 0, 0, 0, 8., 0, 0, 0, 4.;
  xtal::Superlattice superlat(prim->basicstructure->lattice(),
                              xtal::Lattice(S));
  std::shared_ptr<config::Supercell const> tmp_supercell =
      std::make_shared<config::Supercell const>(prim, superlat);
  std::shared_ptr<config::Supercell const> canonical_supercell =
      make_canonical_form(*tmp_supercell);

  Eigen::Matrix3d expected_S;
  expected_S << 4., 0, 0, 0, 4., 0, 0, 0, 8.;
  EXPECT_TRUE(almost_equal(
      canonical_supercell->superlattice.superlattice().lat_column_mat(),
      expected_S));
}

TEST_F(CanonicalFormFCCTest, TestSupercell2) {
  Eigen::Matrix3d S;
  S << 4., 0, 0, 0, 8., 0, 0, 0, 4.;
  xtal::Superlattice superlat(prim->basicstructure->lattice(),
                              xtal::Lattice(S));
  std::shared_ptr<config::Supercell const> tmp_supercell =
      std::make_shared<config::Supercell const>(prim, superlat);
  std::vector<std::shared_ptr<config::Supercell const>> equivalents =
      make_equivalents(*tmp_supercell);

  EXPECT_EQ(equivalents.size(), 3);
}

TEST_F(CanonicalFormFCCTest, Test1) {
  config::Configuration configuration(supercell);
  Eigen::VectorXi &occ = configuration.dof_values.occupation;

  Eigen::VectorXi expected(4);
  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  occ << 0, 0, 1, 0;
  config::Configuration canonical_configuration =
      make_canonical_form(configuration, begin, end);

  expected << 1, 0, 0, 0;
  EXPECT_TRUE(
      almost_equal(canonical_configuration.dof_values.occupation, expected));
}

TEST_F(CanonicalFormFCCTest, Test2) {
  config::Configuration configuration(supercell);
  Eigen::VectorXi &occ = configuration.dof_values.occupation;

  Eigen::VectorXi expected(4);
  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  occ << 0, 0, 1, 0;
  std::vector<config::Configuration> equivalents =
      make_equivalents(configuration, begin, end);

  expected << 0, 0, 0, 1;
  EXPECT_TRUE(almost_equal(equivalents[0].dof_values.occupation, expected));

  expected << 0, 0, 1, 0;
  EXPECT_TRUE(almost_equal(equivalents[1].dof_values.occupation, expected));

  expected << 0, 1, 0, 0;
  EXPECT_TRUE(almost_equal(equivalents[2].dof_values.occupation, expected));

  expected << 1, 0, 0, 0;
  EXPECT_TRUE(almost_equal(equivalents[3].dof_values.occupation, expected));
}

class CanonicalFormFCCTest2 : public testing::Test {
 protected:
  CanonicalFormFCCTest2() {
    std::shared_ptr<config::Prim const> prim =
        config::make_shared_prim(test::FCC_binary_prim());
    Eigen::Matrix3l T;
    T << 4, 0, 0, 0, 1, 0, 0, 0, 1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(CanonicalFormFCCTest2, Test1) {
  config::Configuration configuration(supercell);
  Eigen::VectorXi &occ = configuration.dof_values.occupation;

  Eigen::VectorXi expected(4);
  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  occ << 0, 0, 0, 1;
  std::vector<config::Configuration> equivalents =
      make_equivalents(configuration, begin, end);

  {
    config::SupercellSymOp op = from_canonical(configuration, begin, end);
    EXPECT_EQ(op.supercell_factor_group_index(), 0);
    EXPECT_EQ(op.translation_index(), 3);
  }

  {
    config::SupercellSymOp op = to_canonical(configuration, begin, end);
    EXPECT_EQ(op.supercell_factor_group_index(), 0);
    EXPECT_EQ(op.translation_index(), 1);
  }
}

class CanonicalFormFCCTernaryGLStrainDispTest : public testing::Test {
 protected:
  CanonicalFormFCCTernaryGLStrainDispTest() {
    auto prim =
        config::make_shared_prim(test::FCC_ternary_GLstrain_disp_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(CanonicalFormFCCTernaryGLStrainDispTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues &dof_values = configuration.dof_values;
  dof_values.local_dof_values.at("disp")(2, 2) = 1.0;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  std::vector<config::Configuration> equivalents =
      make_equivalents(configuration, begin, end);

  // +/- (x,y,z), per 4 sites
  EXPECT_EQ(equivalents.size(), 6 * 4);
}

TEST_F(CanonicalFormFCCTernaryGLStrainDispTest, Test2) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues &dof_values = configuration.dof_values;
  dof_values.global_dof_values.at("GLstrain")(2) = 0.01;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  std::vector<config::Configuration> equivalents =
      make_equivalents(configuration, begin, end);

  // Exx,Eyy,Ezz
  EXPECT_EQ(equivalents.size(), 3);
}

TEST_F(CanonicalFormFCCTernaryGLStrainDispTest, Test3) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues &dof_values = configuration.dof_values;
  dof_values.occupation(2) = 1;
  dof_values.local_dof_values.at("disp")(2, 2) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(2) = 0.01;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  config::Configuration canonical_configuration =
      make_canonical_form(configuration, begin, end);

  Eigen::VectorXi expected_occ(4);
  expected_occ << 1, 0, 0, 0;
  EXPECT_TRUE(almost_equal(
      expected_occ, canonical_configuration.dof_values.occupation.transpose()));

  Eigen::MatrixXd expected_disp(3, 4);
  expected_disp << 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  EXPECT_TRUE(almost_equal(
      expected_disp,
      canonical_configuration.dof_values.local_dof_values.at("disp")));

  Eigen::VectorXd expected_GLstrain(6);
  expected_GLstrain << 0.01, 0., 0., 0., 0., 0.;
  EXPECT_TRUE(almost_equal(
      expected_GLstrain,
      canonical_configuration.dof_values.global_dof_values.at("GLstrain")));
}

TEST_F(CanonicalFormFCCTernaryGLStrainDispTest, Test4) {
  // occupation takes precedence over local DoF values in comparisons
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues &dof_values = configuration.dof_values;
  dof_values.occupation(3) = 1;
  dof_values.local_dof_values.at("disp")(2, 2) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(2) = 0.01;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  config::Configuration canonical_configuration =
      make_canonical_form(configuration, begin, end);

  Eigen::VectorXi expected_occ(4);
  expected_occ << 1, 0, 0, 0;
  EXPECT_TRUE(almost_equal(
      expected_occ, canonical_configuration.dof_values.occupation.transpose()));

  Eigen::MatrixXd expected_disp(3, 4);
  expected_disp << 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  EXPECT_TRUE(almost_equal(
      expected_disp,
      canonical_configuration.dof_values.local_dof_values.at("disp")));

  Eigen::VectorXd expected_GLstrain(6);
  expected_GLstrain << 0.01, 0., 0., 0., 0., 0.;
  EXPECT_TRUE(almost_equal(
      expected_GLstrain,
      canonical_configuration.dof_values.global_dof_values.at("GLstrain")));
}
