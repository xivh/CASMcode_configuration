#include "casm/configuration/Prim.hh"

#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

TEST(PrimTest, Test1) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_binary_prim());

  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim->sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim->sym_info.global_dof_symgroup_rep.size(), 0);
}

TEST(PrimTest, Test2) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::SimpleCubic_ising_prim());

  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 96);
  EXPECT_EQ(prim->sym_info.point_group->element.size(), 96);
  EXPECT_EQ(prim->sym_info.unitcellcoord_symgroup_rep.size(), 96);
  EXPECT_EQ(prim->sym_info.occ_symgroup_rep.size(), 96);
  EXPECT_EQ(prim->sym_info.has_aniso_occs, true);
  EXPECT_EQ(prim->sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim->sym_info.global_dof_symgroup_rep.size(), 0);
}

TEST(PrimTest, Test3) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_ternary_GLstrain_prim());

  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim->sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim->sym_info.global_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim->sym_info.global_dof_symgroup_rep.at("GLstrain").size(), 48);
}

TEST(PrimTest, Test4) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::SimpleCubic_disp_prim());

  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim->sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim->sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim->sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim->sym_info.local_dof_symgroup_rep.at("disp").size(), 48);
  EXPECT_EQ(prim->sym_info.global_dof_symgroup_rep.size(), 0);
}
