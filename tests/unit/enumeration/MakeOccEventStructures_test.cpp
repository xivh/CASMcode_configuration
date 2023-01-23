#include "casm/configuration/enumeration/MakeOccEventStructures.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/io/SymInfo_json_io.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

class MakeOccEventStructuresTest : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<occ_events::OccSystem> system;

  MakeOccEventStructuresTest() {
    auto basicstructure = std::make_shared<xtal::BasicStructure const>(
        test::FCC_binary_vacancy_prim());
    prim = std::make_shared<config::Prim>(basicstructure);
    system = std::make_shared<occ_events::OccSystem>(
        prim->basicstructure,
        occ_events::make_chemical_name_list(
            *prim->basicstructure, prim->sym_info.factor_group->element));
  }
};

// basic test: FCC A-B-Va, conventional 4-site cell with 1NN A-Va exchange
TEST_F(MakeOccEventStructuresTest, Test1) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // event: sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, -1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, -1}, "Va", 0),
                      system->make_atom_position({0, 0, 0, 0}, "Va", 0)})});

  // default configuration motif, conventional 4-site FCC
  Eigen::Matrix3d L_motif;
  L_motif.col(0) << 4., 0., 0.;
  L_motif.col(1) << 0., 4., 0.;
  L_motif.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L_motif));
  Configuration motif(motif_supercell);

  // test 5 images
  {
    config::MakeOccEventStructures f(motif, event, system);
    int n = 5;
    for (int i = 0; i < n; ++i) {
      xtal::SimpleStructure structure = f(double(i) / (n - 1));
      // structure.within();
      jsonParser json;
      std::cout << "i: " << i << std::endl;
      std::cout << to_json(structure, json) << std::endl;
      std::cout << std::endl;
      EXPECT_EQ(structure.atom_info.size(), 3);
    }
  }

  // test skip_event_occupants
  {
    bool skip_event_occupants = true;
    config::MakeOccEventStructures f(motif, event, system,
                                     skip_event_occupants);

    xtal::SimpleStructure structure = f(0.0);
    // structure.within();
    jsonParser json;
    std::cout << to_json(structure, json) << std::endl;
    std::cout << std::endl;
    EXPECT_EQ(structure.atom_info.size(), 2);
  }
}
