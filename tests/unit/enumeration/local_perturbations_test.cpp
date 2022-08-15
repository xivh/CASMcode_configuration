#include "casm/configuration/Configuration.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/configuration/enumeration/OccEventInfo.hh"
#include "casm/configuration/enumeration/perturbations.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/io/SymInfo_json_io.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

using namespace CASM;

class FCCBinaryLocalPerturbationsTest : public testing::Test {
 protected:
  FCCBinaryLocalPerturbationsTest() {
    auto basicstructure =
        std::make_shared<xtal::BasicStructure const>(test::FCC_binary_prim());
    prim = std::make_shared<config::Prim>(basicstructure);
    system = std::make_shared<occ_events::OccSystem>(
        prim->basicstructure,
        occ_events::make_chemical_name_list(
            *prim->basicstructure, prim->sym_info.factor_group->element));
  }

  void print_coords(config::Configuration const &config) const {
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

  void print_group(config::SymGroup const &group,
                   xtal::Lattice const &lattice) {
    xtal::SymInfoOptions opt{FRAC};
    for (Index i = 0; i < group.element.size(); ++i) {
      auto const &op = group.element[i];
      Index index = group.head_group_index[i];
      xtal::SymInfo syminfo{op, lattice};
      std::cout << index << ": " << to_brief_unicode(syminfo, opt) << std::endl;
    }
  }

  /// print UnitCellCoord as is (not "within")
  void print_local_orbits_sites(
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits) const {
    jsonParser json;
    json["local_orbits"].put_array();
    for (auto const &local_orbit : local_orbits) {
      jsonParser torbit_json;
      torbit_json.put_array();
      for (auto const &local_cluster : local_orbit) {
        jsonParser tclust_json;
        tclust_json.put_array();
        for (auto const &site : local_cluster) {
          tclust_json.push_back(site);
        }
        torbit_json.push_back(tclust_json);
      }
      json["local_orbits"].push_back(torbit_json);
    }
    std::cout << json << std::endl;
  }

  /// print UnitCellCoord "within" supercell
  void print_local_orbits_sites(
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits,
      xtal::UnitCellCoordIndexConverter const &converter) const {
    jsonParser json;
    json["local_orbits"].put_array();
    for (auto const &local_orbit : local_orbits) {
      jsonParser torbit_json;
      torbit_json.put_array();
      for (auto const &local_cluster : local_orbit) {
        jsonParser tclust_json;
        tclust_json.put_array();
        for (auto const &site : local_cluster) {
          tclust_json.push_back(converter(converter(site)));
        }
        torbit_json.push_back(tclust_json);
      }
      json["local_orbits"].push_back(torbit_json);
    }
    std::cout << json << std::endl;
  }

  void print_local_orbits_site_indices(
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits,
      xtal::UnitCellCoordIndexConverter const &converter) const {
    jsonParser json;
    to_json(make_orbits_as_indices(local_orbits, converter),
            json["local_orbits_as_indices"]);
    std::cout << json << std::endl;
  }

  void print_distinct_cluster_site_indices(
      config::Configuration const &background,
      std::vector<std::set<clust::IntegralCluster>> const &local_orbits,
      config::OccEventSupercellInfo const &event_supercell_info) const {
    auto local_orbits_as_indices = clust::make_orbits_as_indices(
        local_orbits,
        event_supercell_info.supercell->unitcellcoord_index_converter);
    auto distinct_local_cluster_sites =
        config::make_distinct_local_cluster_sites(
            background, event_supercell_info.sites,
            event_supercell_info.occ_init, event_supercell_info.occ_final,
            event_supercell_info.supercellsymop_symgroup_rep,
            local_orbits_as_indices);
    jsonParser json;
    json["distinct_clusters"].put_array();
    for (auto const &local_cluster : distinct_local_cluster_sites) {
      jsonParser tjson;
      tjson.put_array();
      for (auto const &index : local_cluster) {
        tjson.push_back(index);
      }
      json["distinct_clusters"].push_back(tjson);
    }
    std::cout << json << std::endl;
  }

  /// check_perturbations does the following, but with more printing:
  ///
  /// \code
  /// auto event_prim_info = std::make_shared<OccEventPrimInfo>(prim, event);
  /// auto supercell = std::make_shared<Supercell const>(prim,
  /// supercell_lattice); OccEventSupercellInfo
  /// event_supercell_info(event_prim_info, supercell);
  ///
  /// std::set<Configuration> backgrounds =
  ///     event_supercell_info.make_distinct_background_configurations(
  ///         motif);
  /// EXPECT_EQ(backgrounds.size(),
  /// expected_perturbations_by_background.size()); Index i = 0; for (auto const
  /// &background : backgrounds) {
  ///   std::set<Configuration> perturbations =
  ///       event_supercell_info.make_distinct_local_perturbations(
  ///           *it, local_clusters);
  ///   EXPECT_EQ(perturbations.size(),
  ///             expected_perturbations_by_background[i]);
  ///   ++i;
  /// }
  ///
  /// std::set<Configuration> all =
  ///     event_supercell_info.make_all_distinct_local_perturbations(
  ///         motif, local_clusters);
  /// EXPECT_EQ(all.size(), expected_total_perturbations);
  /// \endcode
  ///
  void check_perturbations(
      occ_events::OccEvent const &event,
      std::set<clust::IntegralCluster> const &local_clusters,
      config::Configuration const &motif,
      xtal::Lattice const &supercell_lattice,
      std::vector<Index> const &expected_perturbations_by_background,
      double expected_total_perturbations) {
    using namespace clust;
    using namespace config;
    using namespace occ_events;

    std::cout << "prim factor group:" << std::endl;
    print_group(*prim->sym_info.factor_group, prim->basicstructure->lattice());
    std::cout << std::endl;

    // event info
    auto event_prim_info = std::make_shared<OccEventPrimInfo>(prim, event);
    std::cout << "event group:" << std::endl;
    print_group(*event_prim_info->invariant_group,
                prim->basicstructure->lattice());
    std::cout << std::endl;

    // local-cluster orbits
    auto local_orbits = event_prim_info->make_local_orbits(local_clusters);
    print_local_orbits_sites(local_orbits);

    // supercell where perturbations will be made
    auto supercell = std::make_shared<Supercell const>(prim, supercell_lattice);
    OccEventSupercellInfo event_supercell_info(event_prim_info, supercell);

    std::cout << "supercell lattice: \n"
              << supercell->superlattice.superlattice().lat_column_mat()
              << std::endl;
    std::cout << "supercell factor group:" << std::endl;
    print_group(*supercell->sym_info.factor_group,
                prim->basicstructure->lattice());
    std::cout << std::endl;

    // local-cluster orbits in the supercell
    print_local_orbits_sites(local_orbits,
                             supercell->unitcellcoord_index_converter);
    print_local_orbits_site_indices(local_orbits,
                                    supercell->unitcellcoord_index_converter);

    // --- perturbations by distinct background/event combination: ---

    // generate distinct background/event combinations
    std::set<Configuration> backgrounds =
        event_supercell_info.make_distinct_background_configurations(motif);
    EXPECT_EQ(backgrounds.size(), expected_perturbations_by_background.size());

    // for each, generate local perturbations:
    Index i = 0;
    for (auto it = backgrounds.begin(); it != backgrounds.end(); ++it) {
      // print background and unique local-cluster sites,
      //   taking into account the background configuration
      std::cout << "background." << i << ": "
                << it->dof_values.occupation.transpose() << std::endl;
      print_distinct_cluster_site_indices(*it, local_orbits,
                                          event_supercell_info);

      // generate perturbations
      std::set<Configuration> perturbations =
          event_supercell_info.make_distinct_local_perturbations(
              *it, local_clusters);

      // print perturbations
      std::cout << "distinct local perturbations: " << std::endl;
      for (auto const &c : perturbations) {
        std::cout << c.dof_values.occupation.transpose() << std::endl;
      }
      std::cout << std::endl;

      // check number for this background/event combination
      EXPECT_EQ(perturbations.size(), expected_perturbations_by_background[i]);
      ++i;
    }

    // --- all perturbations ---

    // make all distinct perturbations around an event within given background
    // motif
    std::set<Configuration> all =
        event_supercell_info.make_all_distinct_local_perturbations(
            motif, local_clusters);

    // print perturbations
    std::cout << "all perturbations: " << std::endl;
    for (auto const &c : all) {
      // std::cout << "###" << std::endl;
      std::cout << c.dof_values.occupation.transpose() << std::endl;
      // print_coords(c);
      // std::cout << std::endl;
    }

    // check total number
    EXPECT_EQ(all.size(), expected_total_perturbations);
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<occ_events::OccSystem> system;
};

// basic test
// 12x12x12 default configuration + 1-site local cluster -> 2 local environments
TEST_F(FCCBinaryLocalPerturbationsTest, Test1) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // event: sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "B", 0)})});

  // clusters to perturb:
  std::set<IntegralCluster> local_clusters({IntegralCluster({{0, 1, 0, 0}})});

  // default configuration motif
  Eigen::Matrix3d L_motif;
  L_motif.col(0) << 4., 0., 0.;
  L_motif.col(1) << 0., 4., 0.;
  L_motif.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L_motif));
  Configuration motif(motif_supercell);

  // supercell lattice
  Eigen::Matrix3d L_supercell;
  L_supercell.col(0) << 12., 0., 0.;
  L_supercell.col(1) << 0., 12., 0.;
  L_supercell.col(2) << 0., 0., 12.;
  xtal::Lattice supercell_lattice(L_supercell);

  // expected perturbations:
  std::vector<Index> expected_perturbations_by_background({2});
  double expected_total_perturbations = 2;

  check_perturbations(event, local_clusters, motif, supercell_lattice,
                      expected_perturbations_by_background,
                      expected_total_perturbations);
}

// test symmetry breaking due to supercell
// 12x8x8 default configuration + 1-site local cluster -> 3 local environments
TEST_F(FCCBinaryLocalPerturbationsTest, Test2) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // event: sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "B", 0)})});

  // clusters to perturb:
  std::set<IntegralCluster> local_clusters({IntegralCluster({{0, -1, 1, 1}})});

  // default configuration motif
  Eigen::Matrix3d L_motif;
  L_motif.col(0) << 4., 0., 0.;
  L_motif.col(1) << 0., 4., 0.;
  L_motif.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L_motif));
  Configuration motif(motif_supercell);

  // supercell lattice
  Eigen::Matrix3d L_supercell;
  L_supercell.col(0) << 12., 0., 0.;
  L_supercell.col(1) << 0., 8., 0.;
  L_supercell.col(2) << 0., 0., 8.;
  xtal::Lattice supercell_lattice(L_supercell);

  // expected perturbations:
  std::vector<Index> expected_perturbations_by_background({3});
  double expected_total_perturbations = 3;

  check_perturbations(event, local_clusters, motif, supercell_lattice,
                      expected_perturbations_by_background,
                      expected_total_perturbations);
}

// test symmetry breaking due to background
// 12x12x12 L12 configuration + 1-site local cluster -> 5 local environments
TEST_F(FCCBinaryLocalPerturbationsTest, Test3) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // event: sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "B", 0)})});

  // clusters to perturb:
  std::set<IntegralCluster> local_clusters({IntegralCluster({{0, 1, 0, 0}})});

  // default configuration motif
  Eigen::Matrix3d L_motif;
  L_motif.col(0) << 4., 0., 0.;
  L_motif.col(1) << 0., 4., 0.;
  L_motif.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L_motif));
  Configuration motif(motif_supercell);
  motif.dof_values.occupation << 1, 0, 0, 0;

  // supercell lattice
  Eigen::Matrix3d L_supercell;
  L_supercell.col(0) << 12., 0., 0.;
  L_supercell.col(1) << 0., 12., 0.;
  L_supercell.col(2) << 0., 0., 12.;
  xtal::Lattice supercell_lattice(L_supercell);

  // expected perturbations: 5
  // A-A site: no defects
  // A-A site: A on B
  // A-A site: B on A
  // A-B site: no defects
  // A-B site: B on A
  std::vector<Index> expected_perturbations_by_background({3, 2});
  double expected_total_perturbations = 5;

  check_perturbations(event, local_clusters, motif, supercell_lattice,
                      expected_perturbations_by_background,
                      expected_total_perturbations);
}

// test symmetry breaking due to background
// 12x12x12 L12 configuration + 2-site local cluster -> 7 local environments
TEST_F(FCCBinaryLocalPerturbationsTest, Test4) {
  using namespace clust;
  using namespace config;
  using namespace occ_events;

  // event: sites at origin and xy-face center
  OccEvent event(
      {OccTrajectory({system->make_atom_position({0, 0, 0, 0}, "A", 0),
                      system->make_atom_position({0, 0, 0, 1}, "A", 0)}),
       OccTrajectory({system->make_atom_position({0, 0, 0, 1}, "B", 0),
                      system->make_atom_position({0, 0, 0, 0}, "B", 0)})});

  // clusters to perturb:
  std::set<IntegralCluster> local_clusters(
      {IntegralCluster({{0, 1, 0, 0}, {0, 0, 1, 0}})});

  // default configuration motif
  Eigen::Matrix3d L_motif;
  L_motif.col(0) << 4., 0., 0.;
  L_motif.col(1) << 0., 4., 0.;
  L_motif.col(2) << 0., 0., 4.;
  auto motif_supercell =
      std::make_shared<Supercell const>(prim, xtal::Lattice(L_motif));
  Configuration motif(motif_supercell);
  motif.dof_values.occupation << 1, 0, 0, 0;

  // supercell lattice
  Eigen::Matrix3d L_supercell;
  L_supercell.col(0) << 12., 0., 0.;
  L_supercell.col(1) << 0., 12., 0.;
  L_supercell.col(2) << 0., 0., 12.;
  xtal::Lattice supercell_lattice(L_supercell);

  // expected perturbations: 5
  // A-A site: no defects
  // A-A site: A on B
  // A-A site: B on A
  // A-B site: no defects
  // A-B site: B on A
  std::vector<Index> expected_perturbations_by_background({4, 3});
  double expected_total_perturbations = 7;

  check_perturbations(event, local_clusters, motif, supercell_lattice,
                      expected_perturbations_by_background,
                      expected_total_perturbations);
}
