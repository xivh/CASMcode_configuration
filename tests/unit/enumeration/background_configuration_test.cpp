#include "casm/configuration/Configuration.hh"
#include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

using namespace CASM;

class FCCBinaryBackgroundConfigurationTest : public testing::Test {
 protected:
  FCCBinaryBackgroundConfigurationTest() {
    auto basicstructure =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    prim = std::make_shared<config::Prim>(basicstructure);
    system = std::make_shared<occ_events::OccSystem>(
        prim->basicstructure,
        occ_events::make_chemical_name_list(
            *prim->basicstructure, prim->sym_info.factor_group->element));
  }

  void print_coords(config::Configuration const &config) {
    auto const &basicstructure = *config.supercell->prim->basicstructure;
    auto const &converter = config.supercell->unitcellcoord_index_converter;
    auto const &basis = basicstructure.basis();
    auto const &occupation = config.dof_values.occupation;
    std::cout << "---" << std::endl;
    for (Index l = 0; l < occupation.size(); ++l) {
      auto ucc = converter(l);
      Index b = ucc.sublattice();
      std::string moltype = basis[b].occupant_dof()[occupation[l]].name();
      std::cout << ucc.coordinate(basicstructure).const_cart().transpose()
                << " " << moltype << std::endl;
    }
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<occ_events::OccSystem> system;
};

TEST_F(FCCBinaryBackgroundConfigurationTest, Test1) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "B", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "A", 0)})});
  auto event_prim_info = std::make_shared<OccEventPrimInfo>(prim, event);

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L));
  Configuration motif(motif_supercell);
  motif.dof_values.occupation << 1, 0, 0, 0;

  OccEventSupercellInfo event_supercell_info(event_prim_info, motif_supercell);
  std::set<Configuration> backgrounds =
      event_supercell_info.make_distinct_background_configurations(motif);
  // for (auto const &c : backgrounds) {
  //   std::cout << c.dof_values.occupation.transpose() << std::endl;
  // }
  EXPECT_EQ(backgrounds.size(), 2);
}

TEST_F(FCCBinaryBackgroundConfigurationTest, Test2) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "B", 0)})});
  auto event_prim_info = std::make_shared<OccEventPrimInfo>(prim, event);

  Eigen::Matrix3d L;
  // conventional 4-atom fcc supercell
  L.col(0) << 4., 0., 0.;
  L.col(1) << 0., 4., 0.;
  L.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L));
  Configuration motif(motif_supercell);
  motif.dof_values.occupation << 1, 0, 0, 0;

  // conventional 4-atom fcc supercell
  L.col(0) << 12., 0., 0.;
  L.col(1) << 0., 12., 0.;
  L.col(2) << 0., 0., 12.;
  auto supercell = std::make_shared<Supercell const>(prim, xtal::Lattice(L));
  OccEventSupercellInfo event_supercell_info(event_prim_info, supercell);

  std::set<Configuration> backgrounds =
      event_supercell_info.make_distinct_background_configurations(motif);
  // for (auto const &c : backgrounds) {
  //   std::cout << c.dof_values.occupation.transpose() << std::endl;
  // }
  EXPECT_EQ(backgrounds.size(), 2);
}
