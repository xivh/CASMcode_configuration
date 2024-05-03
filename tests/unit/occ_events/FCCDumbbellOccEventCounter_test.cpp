#include "casm/configuration/group/Group.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SymInfo.hh"
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
class FCCDumbbellOccEventCounterTest : public testing::Test {
 protected:
  std::shared_ptr<xtal::BasicStructure const> prim;
  std::shared_ptr<occ_events::SymGroup const> factor_group;
  std::vector<occ_events::OccEventRep> occevent_symgroup_rep;
  std::shared_ptr<occ_events::OccSystem> system;

  FCCDumbbellOccEventCounterTest() {
    prim = std::make_shared<xtal::BasicStructure const>(this->make_prim());
    factor_group = sym_info::make_factor_group(*prim);
    occevent_symgroup_rep =
        occ_events::make_occevent_symgroup_rep(factor_group->element, *prim);
    system = std::make_shared<occ_events::OccSystem>(
        prim,
        occ_events::make_chemical_name_list(*prim, factor_group->element));
  }

  xtal::BasicStructure make_prim() {
    using namespace CASM;
    using namespace CASM::xtal;

    // lattice vectors as cols
    Eigen::Matrix3d lat;
    lat << 0.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0, 2.0, 0.0;

    BasicStructure struc{Lattice{lat}};
    struc.set_title("FCC_dumbbell");

    double offset = 0.4;
    bool divisible = true;

    Molecule A = Molecule::make_atom("A");
    Molecule A2x("A2_dumbbell",
                 {AtomPosition(Eigen::Vector3d(-offset, 0.0, 0.0), "A"),
                  AtomPosition(Eigen::Vector3d(offset, 0.0, 0.0), "A")},
                 divisible);
    Molecule A2y("A2_dumbbell",
                 {AtomPosition(Eigen::Vector3d(0.0, -offset, 0.0), "A"),
                  AtomPosition(Eigen::Vector3d(0.0, offset, 0.0), "A")},
                 divisible);
    Molecule A2z("A2_dumbbell",
                 {AtomPosition(Eigen::Vector3d(0.0, 0.0, -offset), "A"),
                  AtomPosition(Eigen::Vector3d(0.0, 0.0, offset), "A")},
                 divisible);
    struc.set_unique_names(
        {{"A", "A2_dumbbell.x", "A2_dumbbell.y", "A2_dumbbell.z"}});

    struc.push_back(Site(
        Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), CART),
        std::vector<Molecule>{A, A2x, A2y, A2z}, std::vector<SiteDoFSet>{}));

    return struc;
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

TEST_F(FCCDumbbellOccEventCounterTest, Test1) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0)})});
  // clang-format on

  {
    OccEventCounterParameters params;
    params.required_init_molecule_count = to_VectorXi({2, 0});
    EXPECT_EQ(count_occevents(clusters, params), 0);
  }

  {
    OccEventCounterParameters params;
    params.skip_direct_exchange = false;
    params.required_init_molecule_count = to_VectorXi({2, 0});
    EXPECT_EQ(count_occevents(clusters, params), 1);
  }
}

TEST_F(FCCDumbbellOccEventCounterTest, Test1b) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0)})});
  // clang-format on

  {
    OccEventCounterParameters params;
    params.required_init_orientation_count = to_VectorXi({1, 0, 0, 1});
    params.required_final_orientation_count = to_VectorXi({1, 0, 0, 1});
    EXPECT_EQ(count_occevents(clusters, params), 9);
  }

  {
    OccEventCounterParameters params;
    params.skip_direct_exchange = false;
    params.required_init_orientation_count = to_VectorXi({1, 0, 0, 1});
    params.required_final_orientation_count = to_VectorXi({1, 0, 0, 1});
    EXPECT_EQ(count_occevents(clusters, params), 14);
  }
}

TEST_F(FCCDumbbellOccEventCounterTest, Test1c) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0)})});
  // clang-format on

  {
    OccEventCounterParameters params;
    params.required_init_orientation_count = to_VectorXi({0, 1, 1, 0});
    EXPECT_EQ(count_occevents(clusters, params), 202);
  }

  {
    OccEventCounterParameters params;
    params.skip_direct_exchange = false;
    params.required_init_orientation_count = to_VectorXi({0, 1, 1, 0});
    EXPECT_EQ(count_occevents(clusters, params), 318);
  }
}

TEST_F(FCCDumbbellOccEventCounterTest, Test2) {
  using namespace CASM::occ_events;

  // clang-format off
  std::vector<clust::IntegralCluster> clusters({
      clust::IntegralCluster({
          xtal::UnitCellCoord(0, 0, 0, 0),
          xtal::UnitCellCoord(0, 1, 0, 0)})});
  // clang-format on

  OccEventCounterParameters params;
  params.required_init_orientation_count = to_VectorXi({1, 0, 0, 1});
  params.required_final_orientation_count = to_VectorXi({1, 0, 0, 1});

  std::vector<OccEvent> prototypes = make_prim_periodic_occevent_prototypes(
      system, clusters, occevent_symgroup_rep, params);

  EXPECT_EQ(prototypes.size(), 5);
}
