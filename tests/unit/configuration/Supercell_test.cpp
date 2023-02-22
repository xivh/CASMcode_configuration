#include "casm/configuration/Supercell.hh"

#include "casm/configuration/Prim.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

TEST(SupercellTest, Test1) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  std::shared_ptr<config::Supercell const> supercell =
      std::make_shared<config::Supercell const>(prim, T);
  EXPECT_EQ(supercell->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(supercell->sym_info.translation_permutations->size(), 4);
  EXPECT_EQ(supercell->sym_info.factor_group_permutations.size(), 48);
}
