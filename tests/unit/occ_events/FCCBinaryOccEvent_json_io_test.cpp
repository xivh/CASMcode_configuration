
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/occ_events/io/json/OccEventCounter_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
// #include "casm/casm_io/container/json.hh"
// #include "casm/casm_io/container/stream_io.hh"
// #include "casm/casm_io/json/jsonParser.hh"
// #include "casm/configuration/occ_events/io/json/OccEventCounter_json_io.hh"
// #include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"
// #include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"
// #include
// "casm/configuration/occ_events/io/stream/OccEventCounter_stream_io.hh"
// #include "casm/configuration/occ_events/io/stream/OccEvent_stream_io.hh"
// #include "casm/crystallography/io/SymInfo_json_io.hh"

// params.print_state_info = OccEventCounterStateInfoPrinter(*system);

using namespace CASM;

namespace {
Eigen::VectorXi to_VectorXi(std::vector<int> const &v) {
  Eigen::VectorXi vec(v.size());
  for (Index i = 0; i < v.size(); ++i) {
    vec(i) = v[i];
  }
  return vec;
}

}  // namespace

// FCC binary tests
class FCCBinaryOccEventJsonIOTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<occ_events::SymGroup const> factor_group;
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;
  std::shared_ptr<occ_events::OccSystem> system;

  FCCBinaryOccEventJsonIOTest() {
    prim =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    factor_group = sym_info::make_factor_group(*prim);
    occevent_symgroup_rep =
        occ_events::make_occevent_symgroup_rep(factor_group->element, *prim);
    system = std::make_shared<occ_events::OccSystem>(
        prim,
        occ_events::make_chemical_name_list(*prim, factor_group->element));
  }

  std::string event_json_str() {
    return std::string(R"({
    "trajectories": [
      [
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 0, 0],
          "occupant_index": 1
        },
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 0, 1],
          "occupant_index": 1
        }
      ],
      [
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 0, 1],
          "occupant_index": 1
        },
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 1, 0],
          "occupant_index": 1
        }
      ],
      [
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 1, 0],
          "occupant_index": 1
        },
        {
          "atom_name": "B",
          "atom_position_index": 0,
          "chemical_name": "B",
          "coordinate": [0, 0, 0, 0],
          "occupant_index": 1
        }
      ]
    ]
  })");
  }
};

/// Used to test JSON output
TEST_F(FCCBinaryOccEventJsonIOTest, Test1) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0),
          xtal::UnitCellCoord(0, 0, 1, 0)})});
  // clang-format on

  OccEventCounterParameters params;

  std::vector<std::set<OccEvent>> orbits = make_prim_periodic_occevent_orbits(
      system, clusters, occevent_symgroup_rep, params);

  EXPECT_EQ(orbits.size(), 4);

  OccEventOutputOptions options;
  // options.sym_info_options.print_matrix_tau = true;
  options.include_elements = true;
  options.include_invariant_group = true;
  options.include_equivalence_map = true;

  jsonParser json;
  json.put_array();
  for (auto const &orbit : orbits) {
    jsonParser tjson;
    to_json(orbit, tjson, *system, factor_group, occevent_symgroup_rep,
            options);
    json.push_back(tjson);
  }
  EXPECT_EQ(json.size(), 4);
}

TEST_F(FCCBinaryOccEventJsonIOTest, Test2) {
  using namespace CASM::occ_events;
  jsonParser json = jsonParser::parse(event_json_str());
  InputParser<OccEvent> parser(json, *system);
  EXPECT_TRUE(parser.valid());
  EXPECT_EQ(parser.value->elements().size(), 3);
}
