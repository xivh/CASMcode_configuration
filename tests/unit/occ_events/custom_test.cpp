#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// IO
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"

using namespace CASM;

class CustomEventOrbitTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<occ_events::SymGroup const> factor_group;
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;
  std::unique_ptr<occ_events::OccSystem> system;

  CustomEventOrbitTest() {
    prim = std::make_shared<xtal::BasicStructure const>(
        test::FCC_binary_vacancy_prim());
    factor_group = sym_info::make_factor_group(*prim);
    occevent_symgroup_rep =
        occ_events::make_occevent_symgroup_rep(factor_group->element, *prim);
    system = std::make_unique<occ_events::OccSystem>(
        prim,
        occ_events::make_chemical_name_list(*prim, factor_group->element));
  }
};

TEST_F(CustomEventOrbitTest, Test1) {
  using namespace CASM::occ_events;

  xtal::UnitCellCoord site0(0, 0, 0, 0);
  xtal::UnitCellCoord site1(0, 1, 0, 0);

  OccEvent occ_event(
      {OccTrajectory({system->make_molecule_position(site0, "A"),
                      system->make_molecule_position(site1, "A")}),
       OccTrajectory({system->make_molecule_position(site1, "Va"),
                      system->make_molecule_position(site0, "Va")})});

  std::set<OccEvent> orbit =
      make_prim_periodic_orbit(occ_event, occevent_symgroup_rep);

  std::shared_ptr<SymGroup const> occevent_group = make_occevent_group(
      occ_event, factor_group, prim->lattice().lat_column_mat(),
      occevent_symgroup_rep);

  OccEventOutputOptions options;
  options.include_elements = true;
  options.include_invariant_group = true;
  options.include_equivalence_map = true;

  // --- lambda to print symop descriptions ---
  auto _put_desc = [&](occ_events::SymGroup const &factor_group,
                       jsonParser &j) {
    j.put_array();
    for (auto const &op : factor_group.element) {
      xtal::SymInfo syminfo{op, prim->lattice()};
      j.push_back(to_brief_unicode(syminfo, options.sym_info_options));
    }
    j.set_multiline_array();
  };

  jsonParser json;
  to_json(*prim, json["prim"], FRAC);
  _put_desc(*factor_group, json["factor_group"]);
  to_json(*system, json["event_system"]);
  to_json(orbit, json["orbit"], *system, factor_group, occevent_symgroup_rep,
          options);
  std::cout << json << std::endl;
}
