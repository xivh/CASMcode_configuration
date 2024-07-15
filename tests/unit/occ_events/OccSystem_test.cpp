#include "casm/configuration/occ_events/OccSystem.hh"

#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

namespace {

// invalid molecule naming: equivalent molecules have different names
inline CASM::xtal::BasicStructure FCC_dimer_prim_invalid_naming1() {
  using namespace CASM;
  using namespace CASM::xtal;

  // lattice vectors as cols
  Eigen::Matrix3d lat;
  lat << 0.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 0.0;

  BasicStructure struc{Lattice{lat}};
  struc.set_title("FCC_dimer");

  Molecule A2x("A2.x", {AtomPosition(Eigen::Vector3d(-0.4, 0.0, 0.0), "A"),
                        AtomPosition(Eigen::Vector3d(0.4, 0.0, 0.0), "A")});
  Molecule A2y("A2.y", {AtomPosition(Eigen::Vector3d(0.0, -0.4, 0.0), "A"),
                        AtomPosition(Eigen::Vector3d(0.0, 0.4, 0.0), "A")});
  Molecule A2z("A2.z", {AtomPosition(Eigen::Vector3d(0.0, 0.0, -0.4), "A"),
                        AtomPosition(Eigen::Vector3d(0.0, 0.0, 0.4), "A")});

  struc.push_back(
      Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), CART),
           std::vector<Molecule>{A2x, A2y, A2z}, std::vector<SiteDoFSet>{}));

  return struc;
}

}  // namespace

// valid molecule naming: equivalent (atomic) molecules have the same names
TEST(MoleculeNamingTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());

  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  EXPECT_TRUE(
      occ_events::is_valid_molecule_naming(*prim, factor_group->element));
}

// valid molecule naming: equivalent molecules have the same names
TEST(MoleculeNamingTest, Test2) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());

  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  EXPECT_TRUE(
      occ_events::is_valid_molecule_naming(*prim, factor_group->element));
}

// invalid molecule naming: equivalent molecules have different names
TEST(MoleculeNamingTest, Test3) {
  auto prim = std::make_shared<xtal::BasicStructure const>(
      FCC_dimer_prim_invalid_naming1());

  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  EXPECT_FALSE(
      occ_events::is_valid_molecule_naming(*prim, factor_group->element));
}

TEST(OccSystemTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> chemical_name_list =
      occ_events::make_chemical_name_list(*prim, factor_group->element);

  EXPECT_EQ(chemical_name_list, std::vector<std::string>({"A", "B"}));
}

TEST(OccSystemTest, Test2) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> chemical_name_list =
      occ_events::make_chemical_name_list(*prim, factor_group->element);

  EXPECT_EQ(chemical_name_list, std::vector<std::string>({"A2"}));
}

TEST(OccSystemTest, Test3) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> orientation_name_list =
      occ_events::make_orientation_name_list(*prim);

  EXPECT_EQ(orientation_name_list,
            std::vector<std::string>({"A2.x", "A2.y", "A2.z"}));
}

TEST(MakereservoirPositionTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> orientation_name_list =
      occ_events::make_orientation_name_list(*prim);
  EXPECT_EQ(orientation_name_list,
            std::vector<std::string>({"A2.x", "A2.y", "A2.z"}));

  std::vector<std::string> chemical_name_list =
      occ_events::make_chemical_name_list(*prim, factor_group->element);
  EXPECT_EQ(chemical_name_list, std::vector<std::string>({"A2"}));

  occ_events::OccSystem system(prim, chemical_name_list);
  occ_events::OccPosition occ_position =
      system.make_molecule_in_reservoir_position("A2");

  EXPECT_EQ(occ_position.is_in_reservoir, true);
  EXPECT_EQ(occ_position.is_atom, false);
  EXPECT_EQ(occ_position.integral_site_coordinate,
            xtal::UnitCellCoord(0, 0, 0, 0));
  EXPECT_EQ(occ_position.occupant_index, 0);
  EXPECT_EQ(occ_position.atom_position_index, -1);
}

TEST(MakeMoleculePositionTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> orientation_name_list =
      occ_events::make_orientation_name_list(*prim);
  EXPECT_EQ(orientation_name_list,
            std::vector<std::string>({"A2.x", "A2.y", "A2.z"}));

  std::vector<std::string> chemical_name_list =
      occ_events::make_chemical_name_list(*prim, factor_group->element);
  EXPECT_EQ(chemical_name_list, std::vector<std::string>({"A2"}));

  occ_events::OccSystem system(prim, chemical_name_list);
  occ_events::OccPosition occ_position =
      system.make_molecule_position(xtal::UnitCellCoord(0, 1, 0, 0), "A2.y");

  EXPECT_EQ(occ_position.is_in_reservoir, false);
  EXPECT_EQ(occ_position.is_atom, false);
  EXPECT_EQ(occ_position.integral_site_coordinate,
            xtal::UnitCellCoord(0, 1, 0, 0));
  EXPECT_EQ(occ_position.occupant_index, 1);
  EXPECT_EQ(occ_position.atom_position_index, -1);
}

TEST(MakeAtomicComponentPositionTest, Test1) {
  auto prim =
      std::make_shared<xtal::BasicStructure const>(test::FCC_dimer_prim());
  std::shared_ptr<occ_events::SymGroup const> factor_group =
      sym_info::make_factor_group(*prim);

  std::vector<std::string> orientation_name_list =
      occ_events::make_orientation_name_list(*prim);
  EXPECT_EQ(orientation_name_list,
            std::vector<std::string>({"A2.x", "A2.y", "A2.z"}));

  std::vector<std::string> chemical_name_list =
      occ_events::make_chemical_name_list(*prim, factor_group->element);
  EXPECT_EQ(chemical_name_list, std::vector<std::string>({"A2"}));

  occ_events::OccSystem system(prim, chemical_name_list);
  occ_events::OccPosition occ_position =
      system.make_atom_position(xtal::UnitCellCoord(0, 1, 0, 0), "A2.y", 1);

  EXPECT_EQ(occ_position.is_in_reservoir, false);
  EXPECT_EQ(occ_position.is_atom, true);
  EXPECT_EQ(occ_position.integral_site_coordinate,
            xtal::UnitCellCoord(0, 1, 0, 0));
  EXPECT_EQ(occ_position.occupant_index, 1);
  EXPECT_EQ(occ_position.atom_position_index, 1);
}
