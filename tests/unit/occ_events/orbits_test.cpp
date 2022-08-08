#include "casm/configuration/occ_events/orbits.hh"

#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

// FCC binary tests
class FCCBinaryOccEventOrbitTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<occ_events::SymGroup const> factor_group;
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;
  std::unique_ptr<occ_events::OccSystem> system;

  FCCBinaryOccEventOrbitTest() {
    prim =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    factor_group = sym_info::make_factor_group(*prim);
    occevent_symgroup_rep =
        occ_events::make_occevent_symgroup_rep(factor_group->element, *prim);
    system = std::make_unique<occ_events::OccSystem>(
        prim,
        occ_events::make_chemical_name_list(*prim, factor_group->element));
  }
};

TEST_F(FCCBinaryOccEventOrbitTest, Test1) {
  using namespace CASM::occ_events;

  xtal::UnitCellCoord site0(0, 0, 0, 0);
  xtal::UnitCellCoord site1(0, 1, 0, 0);

  OccEvent occ_event(
      {OccTrajectory({system->make_molecule_position(site0, "B"),
                      system->make_molecule_position(site1, "B")}),
       OccTrajectory({system->make_molecule_position(site1, "A"),
                      system->make_molecule_position(site0, "A")})});

  EXPECT_EQ(occ_event.size(), 2);

  std::set<OccEvent> orbit =
      make_prim_periodic_orbit(occ_event, occevent_symgroup_rep);
  EXPECT_EQ(orbit.size(), 6);

  std::shared_ptr<SymGroup const> occevent_group = make_occevent_group(
      occ_event, factor_group, prim->lattice().lat_column_mat(),
      occevent_symgroup_rep);
  EXPECT_EQ(occevent_group->element.size(), 8);
}

TEST_F(FCCBinaryOccEventOrbitTest, Test2) {
  using namespace CASM::occ_events;

  xtal::UnitCellCoord site0(0, 0, 0, 0);
  xtal::UnitCellCoord site1(0, 1, 1, -1);

  OccEvent occ_event(
      {OccTrajectory({system->make_molecule_position(site0, "B"),
                      system->make_molecule_position(site1, "B")}),
       OccTrajectory({system->make_molecule_position(site1, "A"),
                      system->make_molecule_position(site0, "A")})});

  std::set<OccEvent> orbit =
      make_prim_periodic_orbit(occ_event, occevent_symgroup_rep);
  EXPECT_EQ(orbit.size(), 3);

  std::shared_ptr<SymGroup const> occevent_group = make_occevent_group(
      occ_event, factor_group, prim->lattice().lat_column_mat(),
      occevent_symgroup_rep);
  EXPECT_EQ(occevent_group->element.size(), 16);
}

TEST_F(FCCBinaryOccEventOrbitTest, Test3) {
  using namespace CASM::occ_events;

  xtal::UnitCellCoord site0(0, 0, 0, 0);
  xtal::UnitCellCoord site1(0, 1, 0, 0);
  xtal::UnitCellCoord site2(0, 0, 1, 0);

  OccEvent occ_event(
      {OccTrajectory({system->make_molecule_position(site0, "A"),
                      system->make_molecule_position(site1, "A")}),
       OccTrajectory({system->make_molecule_position(site1, "A"),
                      system->make_molecule_position(site2, "A")}),
       OccTrajectory({system->make_molecule_position(site2, "A"),
                      system->make_molecule_position(site0, "A")})});

  std::set<OccEvent> orbit =
      make_prim_periodic_orbit(occ_event, occevent_symgroup_rep);
  EXPECT_EQ(orbit.size(), 8);

  std::shared_ptr<SymGroup const> occevent_group = make_occevent_group(
      occ_event, factor_group, prim->lattice().lat_column_mat(),
      occevent_symgroup_rep);
  EXPECT_EQ(occevent_group->element.size(), 6);
}

TEST_F(FCCBinaryOccEventOrbitTest, Test4) {
  using namespace CASM::occ_events;

  xtal::UnitCellCoord site0(0, 0, 0, 0);
  xtal::UnitCellCoord site1(0, 1, 0, 0);
  xtal::UnitCellCoord site2(0, 0, 1, 0);

  OccEvent occ_event(
      {OccTrajectory({system->make_molecule_position(site0, "B"),
                      system->make_molecule_position(site1, "B")}),
       OccTrajectory({system->make_molecule_position(site1, "A"),
                      system->make_molecule_position(site2, "A")}),
       OccTrajectory({system->make_molecule_position(site2, "A"),
                      system->make_molecule_position(site0, "A")})});

  std::set<OccEvent> orbit =
      make_prim_periodic_orbit(occ_event, occevent_symgroup_rep);
  EXPECT_EQ(orbit.size(), 24);

  std::shared_ptr<SymGroup const> occevent_group = make_occevent_group(
      occ_event, factor_group, prim->lattice().lat_column_mat(),
      occevent_symgroup_rep);
  EXPECT_EQ(occevent_group->element.size(), 2);
}
