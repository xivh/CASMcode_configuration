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
  Eigen::VectorXi zeros(supercell->unitcellcoord_index_converter.total_sites());
  zeros.setZero();
  EXPECT_TRUE(almost_equal(dof_values.occupation, zeros));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  // ConfigDoFValues copy_apply(SupercellSymOp const &op,
  //                            ConfigDoFValues dof_values);
  dof_values.occupation(0) = 1;
  Index size = dof_values.occupation.size();
  Eigen::VectorXi count = Eigen::VectorXi::Zero(size);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    std::cout << transformed_dof_values.occupation.transpose() << std::endl;
    count += transformed_dof_values.occupation;
  }
  Eigen::VectorXi expected = Eigen::VectorXi::Constant(size, 48);
  EXPECT_TRUE(almost_equal(count, expected));
}

class SupercellSymOpFCCDispTest : public testing::Test {
 protected:
  SupercellSymOpFCCDispTest() {
    std::shared_ptr<config::Prim const> prim =
        std::make_shared<config::Prim const>(test::FCC_binary_disp_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(SupercellSymOpFCCDispTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Eigen::VectorXi zeros(supercell->unitcellcoord_index_converter.total_sites());
  zeros.setZero();
  EXPECT_TRUE(almost_equal(dof_values.occupation, zeros));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  // ConfigDoFValues copy_apply(SupercellSymOp const &op,
  //                            ConfigDoFValues dof_values);
  dof_values.occupation(0) = 1;
  dof_values.local_dof_values.at("disp")(0, 0) = 1.0;
  Eigen::VectorXi count = Eigen::VectorXi::Zero(4);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    std::cout << transformed_dof_values.occupation.transpose() << std::endl;
    std::cout << transformed_dof_values.local_dof_values.at("disp") << std::endl
              << std::endl;
  }
}
