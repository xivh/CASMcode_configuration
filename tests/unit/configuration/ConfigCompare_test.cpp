#include "casm/configuration/ConfigCompare.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class ConfigCompareFCCTest : public testing::Test {
 protected:
  ConfigCompareFCCTest() {
    std::shared_ptr<config::Prim const> prim =
        config::make_shared_prim(test::FCC_binary_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(ConfigCompareFCCTest, Test1) {
  // comparing configurations with the same supercell
  config::Configuration lhs_configuration(supercell);
  Eigen::VectorXi &lhs_occ = lhs_configuration.dof_values.occupation;
  config::Configuration rhs_configuration(supercell);
  Eigen::VectorXi &rhs_occ = rhs_configuration.dof_values.occupation;

  lhs_occ << 0, 0, 1, 0;
  config::ConfigIsEquivalent equal_to_f(lhs_configuration);
  config::ConfigCompare compare_f(equal_to_f);

  rhs_occ << 0, 0, 0, 0;
  EXPECT_FALSE(compare_f(rhs_configuration));
  EXPECT_FALSE(equal_to_f(rhs_configuration));

  rhs_occ << 0, 0, 0, 1;
  EXPECT_FALSE(compare_f(rhs_configuration));
  EXPECT_FALSE(equal_to_f(rhs_configuration));

  rhs_occ << 0, 0, 1, 0;
  EXPECT_FALSE(compare_f(rhs_configuration));
  EXPECT_TRUE(equal_to_f(rhs_configuration));

  rhs_occ << 0, 1, 0, 0;
  EXPECT_TRUE(compare_f(rhs_configuration));
  EXPECT_FALSE(equal_to_f(rhs_configuration));

  rhs_occ << 1, 0, 0, 0;
  EXPECT_TRUE(compare_f(rhs_configuration));
  EXPECT_FALSE(equal_to_f(rhs_configuration));

  rhs_occ << 1, 1, 1, 1;
  EXPECT_TRUE(compare_f(rhs_configuration));
  EXPECT_FALSE(equal_to_f(rhs_configuration));
}

TEST_F(ConfigCompareFCCTest, Test2) {
  // comparing configurations with different supercells
  config::Configuration lhs_configuration(supercell);
  Eigen::VectorXi &lhs_occ = lhs_configuration.dof_values.occupation;

  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  auto rhs_supercell =
      std::make_shared<config::Supercell const>(supercell->prim, T);
  config::Configuration rhs_configuration(rhs_supercell);
  Eigen::VectorXi &rhs_occ = rhs_configuration.dof_values.occupation;

  {
    lhs_occ << 0, 0, 1, 0;
    config::ConfigIsEquivalent equal_to_f(lhs_configuration);
    config::ConfigCompare compare_f(equal_to_f);

    rhs_occ << 0;
    EXPECT_FALSE(compare_f(rhs_configuration));
    EXPECT_FALSE(equal_to_f(rhs_configuration));

    rhs_occ << 1;
    EXPECT_FALSE(compare_f(rhs_configuration));
    EXPECT_FALSE(equal_to_f(rhs_configuration));
  }

  // swap lhs and rhs:
  std::swap(lhs_configuration, rhs_configuration);
  {
    lhs_occ << 0;
    config::ConfigIsEquivalent equal_to_f(lhs_configuration);
    config::ConfigCompare compare_f(equal_to_f);

    rhs_occ << 0, 0, 1, 0;
    EXPECT_TRUE(compare_f(rhs_configuration));
    EXPECT_FALSE(equal_to_f(rhs_configuration));
  }

  {
    lhs_occ << 1;
    config::ConfigIsEquivalent equal_to_f(lhs_configuration);
    config::ConfigCompare compare_f(equal_to_f);

    rhs_occ << 0, 0, 1, 0;
    EXPECT_TRUE(compare_f(rhs_configuration));
    EXPECT_FALSE(equal_to_f(rhs_configuration));
  }
}

TEST_F(ConfigCompareFCCTest, Test3) {
  // comparing configurations with different prim (expect throw)
  config::Configuration lhs_configuration(supercell);
  Eigen::VectorXi &lhs_occ = lhs_configuration.dof_values.occupation;

  std::shared_ptr<config::Prim const> rhs_prim =
      config::make_shared_prim(test::FCC_ternary_prim());
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  auto rhs_supercell = std::make_shared<config::Supercell const>(rhs_prim, T);
  config::Configuration rhs_configuration(rhs_supercell);
  Eigen::VectorXi &rhs_occ = rhs_configuration.dof_values.occupation;

  lhs_occ << 0, 0, 1, 0;
  config::ConfigIsEquivalent equal_to_f(lhs_configuration);
  config::ConfigCompare compare_f(equal_to_f);

  rhs_occ << 0;
  EXPECT_THROW(compare_f(rhs_configuration), std::runtime_error);
  EXPECT_THROW(equal_to_f(rhs_configuration), std::runtime_error);
}

TEST_F(ConfigCompareFCCTest, Test4) {
  config::Configuration configuration(supercell);
  Index size = configuration.dof_values.occupation.size();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  configuration.dof_values.occupation(2) = 1;

  std::set<config::Configuration> configurations;
  for (auto it = begin; it != end; ++it) {
    configurations.emplace(copy_apply(*it, configuration));
  }

  // expect 4 configs, ordered: [0,0,0,1], [0,0,1,0], [0,1,0,0], [1,0,0,0]
  EXPECT_EQ(configurations.size(), 4);
  Index l_expected = 3;
  for (auto const &config : configurations) {
    for (Index l = 0; l < size; ++l) {
      if (l == l_expected) {
        EXPECT_EQ(config.dof_values.occupation(l), 1);
      } else {
        EXPECT_EQ(config.dof_values.occupation(l), 0);
      }
    }
    --l_expected;
  }
}
