#include "casm/configuration/supercell_name.hh"

#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SymTools.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class SupercellNameTest : public testing::Test {
 protected:
  xtal::BasicStructure prim;
  std::vector<xtal::SymOp> factor_group;
  std::vector<xtal::SymOp> crystal_point_group;
  Eigen::Vector3d a;
  Eigen::Vector3d b;
  Eigen::Vector3d c;

  SupercellNameTest()
      : prim(test::FCC_binary_prim()),
        factor_group(xtal::make_factor_group(prim)),
        crystal_point_group(xtal::make_crystal_point_group(
            factor_group, prim.lattice().tol())) {
    std::tie(a, b, c) = prim.lattice().vectors();
  }

  bool is_symmetrically_equivalent(xtal::Lattice superlattice1,
                                   xtal::Lattice superlattice2) {
    auto result = is_equivalent_superlattice(
        superlattice1, superlattice2, crystal_point_group.begin(),
        crystal_point_group.end(), prim.lattice().tol());
    return result.first != crystal_point_group.end();
  }
};

TEST_F(SupercellNameTest, Test1) {
  xtal::Lattice superlattice{a, b, c};
  std::string name = config::make_supercell_name(prim.lattice(), superlattice);
  EXPECT_EQ(name, "SCEL1_1_1_1_0_0_0");

  xtal::Lattice recreated_superlattice =
      config::make_superlattice_from_supercell_name(prim.lattice(), name);
  EXPECT_TRUE(recreated_superlattice == superlattice);
  EXPECT_TRUE(
      is_symmetrically_equivalent(recreated_superlattice, superlattice));
}

TEST_F(SupercellNameTest, Test2) {
  // standard cubic FCC unit cell
  xtal::Lattice superlattice{c + b - a, a - b + c, a + b - c};
  std::string name = config::make_supercell_name(prim.lattice(), superlattice);
  EXPECT_EQ(name, "SCEL4_2_2_1_1_1_0");

  xtal::Lattice recreated_superlattice =
      config::make_superlattice_from_supercell_name(prim.lattice(), name);
  EXPECT_FALSE(recreated_superlattice == superlattice);
  EXPECT_TRUE(
      is_symmetrically_equivalent(recreated_superlattice, superlattice));
}

TEST_F(SupercellNameTest, Test3) {
  // non-standard, but equivalent cubic FCC unit cell
  xtal::Lattice superlattice{c + b - a, a - b + c, (a + b - c) + (c + b - a)};
  std::string name = config::make_supercell_name(prim.lattice(), superlattice);
  EXPECT_EQ(name, "SCEL4_2_2_1_1_1_0");

  xtal::Lattice recreated_superlattice =
      config::make_superlattice_from_supercell_name(prim.lattice(), name);
  EXPECT_FALSE(recreated_superlattice == superlattice);
  EXPECT_TRUE(
      is_symmetrically_equivalent(recreated_superlattice, superlattice));
}
