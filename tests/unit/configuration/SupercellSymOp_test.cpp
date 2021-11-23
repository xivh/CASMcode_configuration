#include "casm/configuration/SupercellSymOp.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class SupercellSymOpFCCTest : public testing::Test {
 protected:
  SupercellSymOpFCCTest() {
    std::shared_ptr<config::Prim const> prim =
        std::make_shared<config::Prim const>(test::FCC_binary_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(SupercellSymOpFCCTest, Test1) {
  EXPECT_EQ(supercell->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(supercell->sym_info.translation_permutations.size(), 4);
  EXPECT_EQ(supercell->sym_info.factor_group_permutations.size(), 48);
}

TEST_F(SupercellSymOpFCCTest, Test2) {
  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  Index count = 0;
  for (auto it = begin; it != end; ++it) {
    ++count;
  }
  EXPECT_EQ(count, 4 * 48);
}

TEST_F(SupercellSymOpFCCTest, Test3) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
  }
  EXPECT_TRUE(almost_equal(occ_count, Eigen::VectorXi::Constant(size, 48)));
}

class SupercellSymOpFCCTernaryGLStrainDispTest : public testing::Test {
 protected:
  SupercellSymOpFCCTernaryGLStrainDispTest() {
    auto prim = std::make_shared<config::Prim const>(
        test::FCC_ternary_GLstrain_disp_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  dof_values.local_dof_values.at("disp")(0, 0) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(0) = 0.01;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);
  Eigen::MatrixXd disp_count = Eigen::MatrixXd::Zero(3, size);
  Eigen::VectorXd GLstrain_count = Eigen::VectorXd::Zero(6);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
    disp_count += transformed_dof_values.local_dof_values.at("disp");
    GLstrain_count += transformed_dof_values.global_dof_values.at("GLstrain");
  }
  EXPECT_TRUE(almost_equal(occ_count, Eigen::VectorXi::Constant(size, 48)));
  EXPECT_TRUE(almost_equal(disp_count, Eigen::MatrixXd::Zero(3, size)));

  // positive stretch gets transformed to positive x,y,z
  Eigen::VectorXd expected(6);
  expected << 0.01, 0.01, 0.01, 0.0, 0.0, 0.0;
  expected *= 48.0 * 4.0 / 3.0;
  EXPECT_TRUE(almost_equal(GLstrain_count, expected));

  // shear gets transformed to +/- xy,yz,xz, averages to 0.0
  dof_values.global_dof_values.at("GLstrain") = Eigen::VectorXd::Zero(6);
  dof_values.global_dof_values.at("GLstrain")(4) = 0.01;
  GLstrain_count = Eigen::VectorXd::Zero(6);
  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    GLstrain_count += transformed_dof_values.global_dof_values.at("GLstrain");
  }
  EXPECT_TRUE(almost_equal(GLstrain_count, Eigen::VectorXd::Zero(6)));
}

class SupercellSymOpSimpleCubicIsingTest : public testing::Test {
 protected:
  SupercellSymOpSimpleCubicIsingTest() {
    auto prim =
        std::make_shared<config::Prim const>(test::SimpleCubic_ising_prim());
    Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(SupercellSymOpSimpleCubicIsingTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
  }
  // Expect 48*8 ops -> 7x down, 1x up; 48*8 ops -> 7x up, 1x down
  EXPECT_TRUE(almost_equal(occ_count,
                           Eigen::VectorXi::Constant(size, 1 * 48 + 7 * 48)));
}
