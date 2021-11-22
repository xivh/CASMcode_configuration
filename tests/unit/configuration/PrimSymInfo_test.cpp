#include "casm/configuration/PrimSymInfo.hh"

#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

TEST(PrimSymInfoTest, Test1) {
  config::PrimSymInfo prim_sym_info(test::FCC_binary_prim());

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.basis_permutation_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 0);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 0);
}
