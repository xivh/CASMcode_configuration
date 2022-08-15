#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
// #include "casm/casm_io/container/json_io.hh"
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
class FCCBinaryOccEventCounterTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<occ_events::SymGroup const> factor_group;
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;
  std::shared_ptr<occ_events::OccSystem> system;

  FCCBinaryOccEventCounterTest() {
    prim =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    factor_group = sym_info::make_factor_group(*prim);
    occevent_symgroup_rep =
        occ_events::make_occevent_symgroup_rep(factor_group->element, *prim);
    system = std::make_shared<occ_events::OccSystem>(
        prim,
        occ_events::make_chemical_name_list(*prim, factor_group->element));
  }

  Index count_occevents(std::vector<clust::IntegralCluster> const &clusters,
                        occ_events::OccEventCounterParameters const &params) {
    occ_events::OccEventCounter counter(system, clusters, params);

    std::vector<occ_events::OccEvent> events;
    while (!counter.is_finished()) {
      events.push_back(counter.value());
      counter.advance();
    }
    return events.size();
  }
};

TEST_F(FCCBinaryOccEventCounterTest, Test1) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0),
          xtal::UnitCellCoord(0, 0, 1, 0)})});
  // clang-format on

  {
    OccEventCounterParameters params;
    EXPECT_EQ(count_occevents(clusters, params), 10);
  }

  {
    OccEventCounterParameters params;
    params.skip_direct_exchange = false;
    EXPECT_EQ(count_occevents(clusters, params), 10);
  }

  {
    OccEventCounterParameters params;
    params.allow_subcluster_events = true;
    EXPECT_EQ(count_occevents(clusters, params), 24);
  }

  {
    OccEventCounterParameters params;
    params.allow_subcluster_events = true;
    params.skip_direct_exchange = false;
    EXPECT_EQ(count_occevents(clusters, params), 36);
  }
}

TEST_F(FCCBinaryOccEventCounterTest, Test2) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0),
          xtal::UnitCellCoord(0, 0, 1, 0)})});
  // clang-format on

  OccEventCounterParameters params;
  // params.save_state_info = true;
  // params.required_init_molecule_count = to_VectorXi({1, 2});

  std::vector<OccEvent> prototypes = make_prim_periodic_occevent_prototypes(
      system, clusters, occevent_symgroup_rep, params);

  EXPECT_EQ(prototypes.size(), 4);
}
