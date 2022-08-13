#include "casm/configuration/copy_configuration.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class CopyConfigurationFCCTest : public testing::Test {
 protected:
  CopyConfigurationFCCTest() {
    prim = config::make_shared_prim(test::FCC_binary_prim());

    xtal::Lattice const &prim_lattice = prim->basicstructure->lattice();
    // prim lattice is:
    // a: 0., 2., 2. (column 0)
    // b: 2., 0., 2. (column 1)
    // c: 2., 2., 0. (column 2)

    Eigen::Matrix3d L;

    // conventional 4-atom fcc supercell
    L.col(0) << 4., 0., 0.;
    L.col(1) << 0., 4., 0.;
    L.col(2) << 0., 0., 4.;
    supercell = std::make_shared<config::Supercell const>(
        prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));

    // 2-atom supercell, both atoms in xy plane
    L.col(0) << 4., 0., 0.;
    L.col(1) << 4., 4., 0.;
    L.col(2) << 0., 2., 2.;
    sub_supercell_xy = std::make_shared<config::Supercell const>(
        prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));

    // 2-atom supercell, both atoms in yz plane
    L.col(0) << 2., 2., 0.;
    L.col(1) << 0., 4., 4.;
    L.col(2) << 0., 0., 4.;
    sub_supercell_yz = std::make_shared<config::Supercell const>(
        prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));
  }

  Index to_site_index(std::shared_ptr<config::Supercell const> const &supercell,
                      xtal::UnitCellCoord unitcellcoord) {
    return supercell->unitcellcoord_index_converter(unitcellcoord);
  }

  int &occ(config::Configuration &config, Index site_index) {
    return config.dof_values.occupation(site_index);
  }

  int &occ(config::Configuration &config, xtal::UnitCellCoord unitcellcoord) {
    return occ(config, to_site_index(config.supercell, unitcellcoord));
  }

  Index total_sites(config::Configuration &config) {
    return config.supercell->unitcellcoord_index_converter.total_sites();
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;
  std::shared_ptr<config::Supercell const> sub_supercell_xy;
  std::shared_ptr<config::Supercell const> sub_supercell_yz;
};

TEST_F(CopyConfigurationFCCTest, Test1) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  {
    // Get the sub-configuration including the z=0 sites
    config::Configuration sub_config =
        copy_configuration(configuration, sub_supercell_xy);
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 0}), 1);
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 1}), 1);
  }

  {
    // Get the sub-configuration including the z=1/2 sites
    config::Configuration sub_config = copy_configuration(
        configuration, sub_supercell_xy, xtal::UnitCell(1, 0, 0));
    EXPECT_EQ(occ(sub_config, {0, 1, 0, 0}), 0);
    EXPECT_EQ(occ(sub_config, {0, 0, 1, 0}), 0);
  }
}

TEST_F(CopyConfigurationFCCTest, IsPrimitiveTest1) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;
  EXPECT_EQ(is_primitive(configuration), false);

  {
    // Get the sub-configuration including the z=0 sites
    config::Configuration sub_config =
        copy_configuration(configuration, sub_supercell_xy);
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 0}), 1);
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 1}), 1);
    EXPECT_EQ(is_primitive(sub_config), false);
  }

  {
    // Get the sub-configuration including the z=1/2 sites
    config::Configuration sub_config = copy_configuration(
        configuration, sub_supercell_xy, xtal::UnitCell(1, 0, 0));
    EXPECT_EQ(occ(sub_config, {0, 1, 0, 0}), 0);
    EXPECT_EQ(occ(sub_config, {0, 0, 1, 0}), 0);
    EXPECT_EQ(is_primitive(sub_config), false);
  }

  {
    // Get the sub-configuration including the x=0 sites
    config::Configuration sub_config =
        copy_configuration(configuration, sub_supercell_yz);
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 0}), 1);
    EXPECT_EQ(occ(sub_config, {0, 1, 0, 0}), 0);
    EXPECT_EQ(is_primitive(sub_config), true);
  }

  {
    // Get the sub-configuration including the x=1/2 sites
    config::Configuration sub_config = copy_configuration(
        configuration, sub_supercell_yz, xtal::UnitCell(0, 0, 1));
    EXPECT_EQ(occ(sub_config, {0, 0, 0, 0}), 1);
    EXPECT_EQ(occ(sub_config, {0, 1, 0, 0}), 0);
    EXPECT_EQ(is_primitive(sub_config), true);
  }
}

TEST_F(CopyConfigurationFCCTest, MakePrimitiveTest1) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  // Make the primitive configuration
  config::Configuration primitive_configuration = make_primitive(configuration);

  // The primitive configuration should have 2 sites,
  // 1 of each occupation value type
  EXPECT_EQ(total_sites(primitive_configuration), 2);
  EXPECT_EQ(occ(primitive_configuration, {0, 0, 0, 0}), 1);
  EXPECT_EQ(occ(primitive_configuration, {0, 1, 0, 0}), 0);
}

TEST_F(CopyConfigurationFCCTest, CopyTransformTest1) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  for (auto it = begin; it != end; ++it) {
    config::Configuration test_a =
        copy_configuration(it.prim_factor_group_index(), it.translation_frac(),
                           configuration, supercell);
    config::Configuration test_b = copy_apply(*it, configuration);
    EXPECT_TRUE(almost_equal(test_a.dof_values.occupation,
                             test_b.dof_values.occupation));
  }
}

TEST_F(CopyConfigurationFCCTest, CopyTransformTest2) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  // This creates a 2-site supercell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration_yz(sub_supercell_yz);
  occ(configuration_yz, {0, 0, 0, 0}) = 1;
  occ(configuration_yz, {0, 1, 0, 0}) = 0;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  // applying the same transformation to either the 4-site or the equivalent
  // 2-site configuration, and copying the result into a 4-site configuration
  // should result in the same occupation
  for (auto it = begin; it != end; ++it) {
    config::Configuration test_a =
        copy_configuration(it.prim_factor_group_index(), it.translation_frac(),
                           configuration_yz, supercell);
    config::Configuration test_b = copy_apply(*it, configuration);
    EXPECT_TRUE(almost_equal(test_a.dof_values.occupation,
                             test_b.dof_values.occupation));
  }
}

TEST_F(CopyConfigurationFCCTest, InCanonicalSupercellTest1) {
  // This creates a 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  config::Configuration config_in_canonical_supercell =
      make_in_canonical_supercell(configuration);

  EXPECT_EQ(*config_in_canonical_supercell.supercell, *supercell);

  Eigen::VectorXi expected_occ(4);
  expected_occ << 1, 1, 0, 0;
  EXPECT_TRUE(almost_equal(config_in_canonical_supercell.dof_values.occupation,
                           expected_occ));
}

TEST_F(CopyConfigurationFCCTest, InCanonicalSupercellTest2) {
  // start with non-canonical conventional 4-atom fcc supercell
  Eigen::Matrix3d L;
  L.col(0) << 0., 4., 0.;
  L.col(1) << -4., 0., 0.;
  L.col(2) << 0., 0., 4.;
  xtal::Lattice const &prim_lattice = prim->basicstructure->lattice();
  auto non_canonical_supercell = std::make_shared<config::Supercell const>(
      prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));

  // This creates a non-canonical 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(non_canonical_supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  config::Configuration config_in_canonical_supercell =
      make_in_canonical_supercell(configuration);

  EXPECT_EQ(*config_in_canonical_supercell.supercell, *supercell);

  Eigen::VectorXi expected_occ(4);
  expected_occ << 1, 1, 0, 0;
  EXPECT_TRUE(almost_equal(config_in_canonical_supercell.dof_values.occupation,
                           expected_occ));
}

TEST_F(CopyConfigurationFCCTest, InCanonicalSupercellTest3) {
  // start with non-canonical conventional 4-atom fcc supercell
  Eigen::Matrix3d L;
  L.col(0) << 0., 4., 0.;
  L.col(1) << -4., 0., 0.;
  L.col(2) << 0., 0., 4.;
  xtal::Lattice const &prim_lattice = prim->basicstructure->lattice();
  auto non_canonical_supercell = std::make_shared<config::Supercell const>(
      prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));

  // This creates a non-canonical 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration configuration(non_canonical_supercell);
  occ(configuration, {0, 0, 0, 0}) = 1;
  occ(configuration, {0, 0, 0, 1}) = 1;

  config::Configuration canonical_config =
      make_in_canonical_supercell(make_primitive(configuration));

  EXPECT_TRUE(is_canonical(*canonical_config.supercell));

  xtal::Lattice const &superlat =
      canonical_config.supercell->superlattice.superlattice();
  Eigen::Matrix3d expected_L;
  expected_L.col(0) << 2., 2., 0.;
  expected_L.col(1) << 2., -2., 0.;
  expected_L.col(2) << 0., 0., -4.;
  EXPECT_TRUE(almost_equal(superlat.lat_column_mat(), expected_L));

  Eigen::VectorXi expected_occ(2);
  expected_occ << 1, 0;
  EXPECT_TRUE(
      almost_equal(canonical_config.dof_values.occupation, expected_occ));
}

// --- CopyConfigurationFCCTernaryGLStrainDispTest ---

class CopyConfigurationFCCTernaryGLStrainDispTest : public testing::Test {
 protected:
  CopyConfigurationFCCTernaryGLStrainDispTest() {
    prim = config::make_shared_prim(test::FCC_ternary_GLstrain_disp_prim());
  }

  Index to_site_index(std::shared_ptr<config::Supercell const> const &supercell,
                      xtal::UnitCellCoord unitcellcoord) {
    return supercell->unitcellcoord_index_converter(unitcellcoord);
  }

  int &occ(config::Configuration &config, Index site_index) {
    return config.dof_values.occupation(site_index);
  }

  int &occ(config::Configuration &config, xtal::UnitCellCoord unitcellcoord) {
    return occ(config, to_site_index(config.supercell, unitcellcoord));
  }

  Eigen::DenseBase<Eigen::Matrix<double, -1, -1, 0>>::ColXpr disp(
      config::Configuration &config, Index site_index) {
    return config.dof_values.local_dof_values.at("disp").col(site_index);
  }

  Eigen::DenseBase<Eigen::Matrix<double, -1, -1, 0>>::ColXpr disp(
      config::Configuration &config, xtal::UnitCellCoord unitcellcoord) {
    return disp(config, to_site_index(config.supercell, unitcellcoord));
  }

  Eigen::VectorXd &strain(config::Configuration &config) {
    return config.dof_values.global_dof_values.at("GLstrain");
  }

  void print_one(config::Configuration const &tconfig) {
    std::cout << "occ: " << tconfig.dof_values.occupation.transpose()
              << std::endl;
    std::cout << "disp: \n"
              << tconfig.dof_values.local_dof_values.at("disp") << std::endl;
    std::cout << "strain: \n"
              << tconfig.dof_values.global_dof_values.at("GLstrain").transpose()
              << std::endl;
  }

  void print(config::Configuration const &init_config) {
    std::cout << "----------0" << std::endl;
    auto begin = config::SupercellSymOp::begin(init_config.supercell);
    auto end = config::SupercellSymOp::end(init_config.supercell);
    for (auto it = begin; it != end; ++it) {
      config::Configuration tconfig = copy_apply(*it, init_config);
      print_one(tconfig);
    }
  }

  std::shared_ptr<config::Prim const> prim;
};

TEST_F(CopyConfigurationFCCTernaryGLStrainDispTest, Test1) {
  // start with non-canonical conventional 4-atom fcc supercell
  Eigen::Matrix3d L;
  L.col(0) << 0., 4., 0.;
  L.col(1) << -4., 0., 0.;
  L.col(2) << 0., 0., 4.;
  xtal::Lattice const &prim_lattice = prim->basicstructure->lattice();
  auto non_canonical_supercell = std::make_shared<config::Supercell const>(
      prim, xtal::Superlattice(prim_lattice, xtal::Lattice(L)));

  // This creates a non-canonical 4-site conventional FCC cell,
  // where z=0 has occ=1, z=1/2 has occ=0
  config::Configuration init_config(non_canonical_supercell);
  occ(init_config, {0, 0, 0, 0}) = 1;
  occ(init_config, {0, 0, 0, 1}) = 1;
  disp(init_config, {0, 0, 0, 0})(2) = 0.1;
  disp(init_config, {0, 0, 0, 1})(2) = 0.1;
  strain(init_config)(2) = 0.1;

  {
    auto begin = config::SupercellSymOp::begin(init_config.supercell);
    auto end = config::SupercellSymOp::end(init_config.supercell);
    config::Configuration tconfig =
        make_canonical_form(init_config, begin, end);

    EXPECT_FALSE(is_canonical(*tconfig.supercell));
    EXPECT_TRUE(is_canonical(tconfig, begin, end));
  }

  {
    config::Configuration tconfig = make_in_canonical_supercell(init_config);

    EXPECT_TRUE(is_canonical(*tconfig.supercell));

    auto begin = config::SupercellSymOp::begin(tconfig.supercell);
    auto end = config::SupercellSymOp::end(tconfig.supercell);
    EXPECT_TRUE(is_canonical(tconfig, begin, end));
  }

  {
    config::Configuration tconfig =
        make_in_canonical_supercell(make_primitive(init_config));

    EXPECT_TRUE(is_canonical(*tconfig.supercell));

    auto begin = config::SupercellSymOp::begin(tconfig.supercell);
    auto end = config::SupercellSymOp::end(tconfig.supercell);
    EXPECT_TRUE(is_canonical(tconfig, begin, end));
  }
}
