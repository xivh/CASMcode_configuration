#include "casm/configuration/make_simple_structure.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/StrainConverter.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class MakeSimpleStructureTest : public testing::Test {
 protected:
  MakeSimpleStructureTest() {
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
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(MakeSimpleStructureTest, Test1) {
  config::Configuration configuration(supercell);
  auto &dof_values = configuration.dof_values;
  auto &occupation = dof_values.occupation;

  // occ DoF
  occupation(0) = 1;
  occupation(1) = 1;

  config::ToAtomicStructure f;
  xtal::SimpleStructure structure = f(configuration);

  jsonParser json;
  to_json(structure, json);
  //  std::cout << json << std::endl;

  EXPECT_EQ(structure.atom_info.names.size(), 4);
  EXPECT_EQ(structure.atom_info.coords.rows(), 3);
  EXPECT_EQ(structure.atom_info.coords.cols(), 4);
}

TEST_F(MakeSimpleStructureTest, Test2) {
  config::Configuration configuration(supercell);
  auto &dof_values = configuration.dof_values;
  auto &occupation = dof_values.occupation;

  // occ DoF
  occupation(0) = 1;
  occupation(1) = 1;

  std::map<std::string, Eigen::MatrixXd> local_properties;
  std::map<std::string, Eigen::VectorXd> global_properties;

  // disp property
  Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(3, 4);
  disp.col(0) << 0.01, 0.0, -0.01;
  local_properties.emplace("disp", disp);

  // strain property
  DoFKey strain_key = "Ustrain";
  xtal::StrainConverter DoFstrain_converter(strain_key,
                                            Eigen::MatrixXd::Identity(6, 6));

  Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
  F(2, 2) = 1.2;
  Eigen::VectorXd Ustrain = DoFstrain_converter.from_F(F);
  global_properties.emplace(strain_key, Ustrain);

  // make SimpleStructure
  config::ToAtomicStructure f;
  xtal::SimpleStructure structure =
      f(configuration, local_properties, global_properties);

  jsonParser json;
  to_json(structure, json);
  //  std::cout << json << std::endl;

  // check dimensions
  EXPECT_EQ(structure.atom_info.names.size(), 4);
  EXPECT_EQ(structure.atom_info.coords.rows(), 3);
  EXPECT_EQ(structure.atom_info.coords.cols(), 4);

  // check properties size
  EXPECT_EQ(structure.atom_info.properties.size(), 1);
  EXPECT_EQ(structure.properties.size(), 1);

  // check lattice vectors
  Eigen::Matrix3d const &L_structure = structure.lat_column_mat;
  EXPECT_TRUE(almost_equal(L_structure(0, 0), 4.0));
  EXPECT_TRUE(almost_equal(L_structure(1, 1), 4.0));
  EXPECT_TRUE(almost_equal(L_structure(2, 2), 4.8));

  // check atom coords
  Eigen::MatrixXd const &coords_structure = structure.atom_info.coords;
  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(3, 4);
  expected.col(0) << 0.01, 0.0, -0.012;
  expected.col(1) << 2.0, 0.0, 2.4;
  expected.col(2) << 0.0, 2.0, 2.4;
  expected.col(3) << 2.0, 2.0, 0.0;
  EXPECT_TRUE(almost_equal(coords_structure, expected));

  // check strain
  EXPECT_TRUE(structure.properties.count("Ustrain"));
  Eigen::VectorXd const &Ustrain_structure = structure.properties.at("Ustrain");
  EXPECT_TRUE(almost_equal(Ustrain_structure, Ustrain));

  // check disp
  EXPECT_TRUE(structure.atom_info.properties.count("disp"));
  Eigen::MatrixXd const &disp_structure =
      structure.atom_info.properties.at("disp");
  EXPECT_TRUE(almost_equal(disp_structure, disp));
}

class MakeSimpleStructureTestStrainDisp : public testing::Test {
 protected:
  MakeSimpleStructureTestStrainDisp() {
    prim = config::make_shared_prim(test::FCC_ternary_GLstrain_disp_prim());

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
  }

  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;
};

TEST_F(MakeSimpleStructureTestStrainDisp, Test1) {
  config::Configuration configuration(supercell);
  auto &dof_values = configuration.dof_values;
  auto &occupation = dof_values.occupation;
  auto &local_dof_values = dof_values.local_dof_values;
  auto &global_dof_values = dof_values.global_dof_values;

  // occ DoF
  occupation(0) = 1;
  occupation(1) = 1;

  // disp DoF
  Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(3, 4);
  disp.col(0) << 0.01, 0.0, -0.01;
  local_dof_values.at("disp") = disp;

  // strain DoF
  DoFKey strain_dof_key = "GLstrain";

  xtal::StrainConverter DoFstrain_converter(
      strain_dof_key, prim->global_dof_info.at(strain_dof_key).basis());

  Eigen::Matrix3d F = Eigen::Matrix3d::Identity();
  F(2, 2) = 1.2;
  global_dof_values.at(strain_dof_key) = DoFstrain_converter.from_F(F);

  // make SimpleStructure
  xtal::SimpleStructure structure = make_simple_structure(configuration);

  jsonParser json;
  to_json(structure, json);
  //  std::cout << json << std::endl;

  // check dimensions
  EXPECT_EQ(structure.atom_info.names.size(), 4);
  EXPECT_EQ(structure.atom_info.coords.rows(), 3);
  EXPECT_EQ(structure.atom_info.coords.cols(), 4);

  // check properties size
  EXPECT_EQ(structure.atom_info.properties.size(), 0);
  EXPECT_EQ(structure.properties.size(), 0);

  // check lattice vectors
  Eigen::Matrix3d const &L_structure = structure.lat_column_mat;
  EXPECT_TRUE(almost_equal(L_structure(0, 0), 4.0));
  EXPECT_TRUE(almost_equal(L_structure(1, 1), 4.0));
  EXPECT_TRUE(almost_equal(L_structure(2, 2), 4.8));

  // check atom coords
  Eigen::MatrixXd const &coords_structure = structure.atom_info.coords;
  Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(3, 4);
  expected.col(0) << 0.01, 0.0, -0.012;
  expected.col(1) << 2.0, 0.0, 2.4;
  expected.col(2) << 0.0, 2.0, 2.4;
  expected.col(3) << 2.0, 2.0, 0.0;
  EXPECT_TRUE(almost_equal(coords_structure, expected));

  // check strain
  EXPECT_FALSE(structure.properties.count("Ustrain"));

  // check disp
  EXPECT_FALSE(structure.atom_info.properties.count("disp"));
}
