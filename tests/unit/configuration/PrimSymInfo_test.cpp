#include "casm/configuration/PrimSymInfo.hh"

#include "gtest/gtest.h"
#include "teststructures.hh"

// // -- io to inspect atom_position_symgroup_rep --
// #include "casm/casm_io/container/stream_io.hh"
// #include "casm/crystallography/SymInfo.hh"
// #include "casm/crystallography/io/SymInfo_stream_io.hh"

using namespace CASM;

TEST(PrimSymInfoTest, Test1) {
  config::PrimSymInfo prim_sym_info(test::FCC_binary_prim());

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.atom_position_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 0);
}

TEST(PrimSymInfoTest, Test2) {
  config::PrimSymInfo prim_sym_info(test::SimpleCubic_ising_prim());

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 96);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 96);
  EXPECT_EQ(prim_sym_info.unitcellcoord_symgroup_rep.size(), 96);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 96);
  EXPECT_EQ(prim_sym_info.atom_position_symgroup_rep.size(), 96);
  EXPECT_EQ(prim_sym_info.has_aniso_occs, true);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 0);
}

TEST(PrimSymInfoTest, Test3) {
  config::PrimSymInfo prim_sym_info(test::FCC_ternary_GLstrain_prim());

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.atom_position_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.at("GLstrain").size(), 48);
}

TEST(PrimSymInfoTest, Test4) {
  config::PrimSymInfo prim_sym_info(test::SimpleCubic_disp_prim());

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.atom_position_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.has_aniso_occs, false);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.at("disp").size(), 48);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 0);
}

TEST(PrimSymInfoTest, Test5) {
  auto prim = test::FCC_dimer_prim();
  config::PrimSymInfo prim_sym_info(prim);

  EXPECT_EQ(prim_sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.point_group->element.size(), 48);
  EXPECT_EQ(prim_sym_info.unitcellcoord_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.occ_symgroup_rep.size(), 48);
  EXPECT_EQ(prim_sym_info.atom_position_symgroup_rep.size(), 48);

  // xtal::SymInfoOptions opt{CART};
  for (Index s = 0; s < prim_sym_info.factor_group->element.size(); ++s) {
    auto const &atom_position_rep = prim_sym_info.atom_position_symgroup_rep[s];
    EXPECT_EQ(atom_position_rep.size(), 1);  // 1 sublattice
    EXPECT_EQ(atom_position_rep[0].size(),
              3);  // 3 occupants allowed on sublattice 0
    EXPECT_EQ(atom_position_rep[0][0].size(), 2);  // 2 atom positions in A2x
    EXPECT_EQ(atom_position_rep[0][1].size(), 2);  // 2 atom positions in A2x
    EXPECT_EQ(atom_position_rep[0][2].size(), 2);  // 2 atom positions in A2x
    // std::cout << "---" << std::endl;
    // xtal::SymInfo syminfo{prim_sym_info.factor_group->element[s],
    // prim.lattice()}; std::cout << to_brief_unicode(syminfo, opt) <<
    // std::endl; std::cout << atom_position_rep[0][0] << std::endl; std::cout
    // << atom_position_rep[0][1] << std::endl; std::cout <<
    // atom_position_rep[0][2] << std::endl;
  }

  EXPECT_EQ(prim_sym_info.has_aniso_occs, true);
  EXPECT_EQ(prim_sym_info.local_dof_symgroup_rep.size(), 1);
  EXPECT_EQ(prim_sym_info.global_dof_symgroup_rep.size(), 0);
}
