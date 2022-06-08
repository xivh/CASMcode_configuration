#include "casm/configuration/DoFSpace.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;
using namespace test;

namespace {

Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

}  // namespace

class DoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;  // conventional unit cell

  DoFSpaceTest()
      : prim(std::make_shared<config::Prim const>(
            test::FCC_ternary_GLstrain_disp_prim())),
        supercell(std::make_shared<config::Supercell const>(
            prim, _fcc_conventional_transf_mat())) {}
};

TEST_F(DoFSpaceTest, ConstructorTest1_GLstrain) {
  // Construct the GLstrain DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "GLstrain";
  config::DoFSpace dof_space = config::make_dof_space(dof_key, prim);
  EXPECT_EQ(dof_space.dim, 6);
  EXPECT_EQ(dof_space.dim, dof_space.basis.rows());
}

TEST_F(DoFSpaceTest, ConstructorTest2_disp) {
  // Construct the disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space = config::make_dof_space(dof_key, *supercell);
  EXPECT_EQ(dof_space.dim, 4 * 3);
  EXPECT_EQ(dof_space.dim, dof_space.basis.rows());
}

TEST_F(DoFSpaceTest, ExcludeHomogeneousModeSpace) {
  // Construct the disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = config::make_dof_space(dof_key, *supercell);

  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  EXPECT_EQ(homogeneous_mode_space.rows(), dof_space_0.dim);
  EXPECT_EQ(homogeneous_mode_space.cols(), 3);

  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  EXPECT_EQ(dof_space_1.dim, dof_space_0.dim);
  EXPECT_EQ(dof_space_1.subspace_dim, dof_space_0.subspace_dim - 3);
}

/// Tests on a structure with a restricted local basis (2d displacements), but
/// the basis is the same on all sites
class RestrictedLocalDoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest()
      : prim(std::make_shared<config::Prim const>(make_prim())),
        supercell(std::make_shared<config::Supercell const>(
            prim, _fcc_conventional_transf_mat())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest::make_prim() {
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1", "d2"},            // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{0.000000000000, 1.754750223661, 1.754750223661},
              Eigen::Vector3d{1.754750223661, 0.000000000000, 1.754750223661},
              Eigen::Vector3d{1.754750223661, 1.754750223661, 0.000000000000}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xy}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest, FactorGroupSize) {
  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 16);
}

TEST_F(RestrictedLocalDoFSpaceTest, ExcludeHomogeneousModeSpace) {
  // Construct the restricted disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = config::make_dof_space(dof_key, *supercell);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << dof_space_0.basis << std::endl;
  EXPECT_EQ(dof_space_0.dim, 8);
  EXPECT_EQ(dof_space_0.subspace_dim, 8);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 8);
  EXPECT_EQ(homogeneous_mode_space.cols(), 2);

  // check exclude homogeneous mode space
  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << dof_space_1.basis << std::endl;
  EXPECT_EQ(dof_space_1.dim, 8);
  EXPECT_EQ(dof_space_1.subspace_dim, 6);
}

/// Tests on a structure with a restricted local basis (2d displacements), and
/// the basis is different on from site to site, rigid translations are only
/// possible along z
class RestrictedLocalDoFSpaceTest2 : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest2()
      : prim(std::make_shared<config::Prim const>(make_prim())),
        supercell(std::make_shared<config::Supercell const>(
            prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest2::make_prim() {
  // BCC base structure,
  // - with corner atoms {A, B} allowed to displace in xz
  // - with body-centered atoms {A, B} allowed to displace in yz,
  // -> test case where all atoms may displace,
  //    but rigid translations are only allowed in 1d (z)
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_yz}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest2, FactorGroupSize) {
  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 16);
}

TEST_F(RestrictedLocalDoFSpaceTest2, ExcludeHomogeneousModeSpace) {
  // Construct the restricted disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = config::make_dof_space(dof_key, *supercell);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << dof_space_0.basis << std::endl;
  EXPECT_EQ(dof_space_0.dim, 4);
  EXPECT_EQ(dof_space_0.subspace_dim, 4);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 4);
  EXPECT_EQ(homogeneous_mode_space.cols(), 1);

  // check exclude homogeneous mode space
  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << dof_space_1.basis << std::endl;
  EXPECT_EQ(dof_space_1.dim, 4);
  EXPECT_EQ(dof_space_1.subspace_dim, 3);
}

/// Tests on a structure with a restricted local basis (2d displacements), and
/// the basis is different on from site to site, rigid translations are only
/// possible along z, occupation is different on the two sublattices
class RestrictedLocalDoFSpaceTest3 : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  RestrictedLocalDoFSpaceTest3()
      : prim(std::make_shared<config::Prim const>(make_prim())),
        supercell(std::make_shared<config::Supercell const>(
            prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure RestrictedLocalDoFSpaceTest3::make_prim() {
  // BCC base structure,
  // - with corner atoms {A, B} allowed to displace in xz
  // - with body-centered atoms {C, D} allowed to displace in yz,
  // -> test case where all atoms may displace,
  //    but rigid translations are only allowed in 1d (z),
  //    and the asymmetric unit > 1
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  Molecule C = Molecule::make_atom("C");
  Molecule D = Molecule::make_atom("D");

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {C, D}, {disp_yz}}});
  return struc;
}

TEST_F(RestrictedLocalDoFSpaceTest3, FactorGroupSize) {
  // SymInfoOptions opt{CART};
  // brief_description(log(), factor_group, prim->lattice(), opt);

  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 8);
}

TEST_F(RestrictedLocalDoFSpaceTest3, ExcludeHomogeneousModeSpace) {
  // Construct the restricted disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = config::make_dof_space(dof_key, *supercell);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << dof_space_0.basis << std::endl;
  EXPECT_EQ(dof_space_0.dim, 4);
  EXPECT_EQ(dof_space_0.subspace_dim, 4);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 4);
  EXPECT_EQ(homogeneous_mode_space.cols(), 1);

  // check exclude homogeneous mode space
  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << dof_space_1.basis << std::endl;
  EXPECT_EQ(dof_space_1.dim, 4);
  EXPECT_EQ(dof_space_1.subspace_dim, 3);
}

/// Tests on a structure with a restricted local basis (2d displacements), and
/// the basis is different on from site to site such that no rigid translations
/// are possible
class VariableLocalDoFSpaceTest1 : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;
  std::shared_ptr<config::Supercell const> supercell;  // conventional unit cell

  static xtal::BasicStructure make_prim();

  VariableLocalDoFSpaceTest1()
      : prim(std::make_shared<config::Prim const>(make_prim())),
        supercell(std::make_shared<config::Supercell const>(
            prim, Eigen::Matrix3l::Identity())) {}
};

xtal::BasicStructure VariableLocalDoFSpaceTest1::make_prim() {
  // FCC base structure,
  // - with corner atoms {A, B} allowed to displace in 3d
  // - with face atoms {C, D} allowed to displace in the face plane,
  // -> test case where all atoms may displace, but along different dimensions
  //    so no rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  Molecule C = Molecule::make_atom("C");
  Molecule D = Molecule::make_atom("D");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),       // AnisoVal type
      {"d0x", "d0y", "d0z"},        // axes names
      Eigen::Matrix3d::Identity(),  // basis
      {}                            // excluded_occs
  };

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1x", "d1y"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d2y", "d2z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d3x", "d3z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.0, lat, FRAC}, {C, D}, {disp_xy}},
       Site{Coordinate{0.0, 0.5, 0.5, lat, FRAC}, {C, D}, {disp_yz}},
       Site{Coordinate{0.5, 0.0, 0.5, lat, FRAC}, {C, D}, {disp_xz}}});
  return struc;
}

TEST_F(VariableLocalDoFSpaceTest1, FactorGroupSize) {
  EXPECT_EQ(prim->sym_info.factor_group->element.size(), 48);
}

TEST_F(VariableLocalDoFSpaceTest1, ExcludeHomogeneousModeSpace) {
  // In this structure, all sites allow displacements, but no rigid
  // translations are possible

  // Construct the restricted disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = config::make_dof_space(dof_key, *supercell);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << dof_space_0.basis << std::endl;
  EXPECT_EQ(dof_space_0.dim, 9);
  EXPECT_EQ(dof_space_0.subspace_dim, 9);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << homogeneous_mode_space << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), 9);
  EXPECT_EQ(homogeneous_mode_space.cols(), 0);

  // check exclude homogeneous mode space
  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << dof_space_1.basis << std::endl;
  EXPECT_EQ(dof_space_1.dim, 9);
  EXPECT_EQ(dof_space_1.subspace_dim, 9);
}

/// This test class uses the pattern of the previous tests to allow for
/// customization to tests various structures that are found to be problematic
class DebugLocalDoFSpaceTest : public testing::Test {
 protected:
  std::shared_ptr<config::Prim const> prim;            // must make in test
  std::shared_ptr<config::Supercell const> supercell;  // must make in test

  DebugLocalDoFSpaceTest() {}

  void check_FactorGroupSize(Index factor_group_size);

  void check_ExcludeHomogeneousModeSpace(
      std::pair<Index, Index> initial_dof_space_shape,
      std::pair<Index, Index> homogeneous_mode_space_shape,
      std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes,
      std::pair<Index, Index> symmetry_adapted_dof_space_shape);
};

void DebugLocalDoFSpaceTest::check_FactorGroupSize(Index factor_group_size) {
  EXPECT_EQ(prim->sym_info.factor_group->element.size(), factor_group_size);
}

void DebugLocalDoFSpaceTest::check_ExcludeHomogeneousModeSpace(
    std::pair<Index, Index> initial_dof_space_shape,
    std::pair<Index, Index> homogeneous_mode_space_shape,
    std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes,
    std::pair<Index, Index> symmetry_adapted_dof_space_shape) {
  // In this structure, all sites allow displacements, but no rigid
  // translations are possible

  // print_local_dof_symreps(supercell->sym_info());

  // Construct the restricted disp DoF space.
  config::Configuration config{supercell};
  DoFKey dof_key = "disp";
  config::DoFSpace dof_space_0 = make_dof_space(dof_key, *supercell);
  // std::cout << "including homogeneous_mode_space: \n"
  //           << pretty(dof_space_0.basis) << std::endl;
  EXPECT_EQ(dof_space_0.dim, initial_dof_space_shape.first);
  EXPECT_EQ(dof_space_0.subspace_dim, initial_dof_space_shape.second);

  // check make homogeneous mode space
  Eigen::MatrixXd homogeneous_mode_space =
      config::make_homogeneous_mode_space(dof_space_0);
  // std::cout << "homogeneous_mode_space: \n"
  //           << pretty(homogeneous_mode_space) << std::endl;
  EXPECT_EQ(homogeneous_mode_space.rows(), homogeneous_mode_space_shape.first);
  EXPECT_EQ(homogeneous_mode_space.cols(), homogeneous_mode_space_shape.second);

  // check exclude homogeneous mode space
  config::DoFSpace dof_space_1 =
      config::exclude_homogeneous_mode_space(dof_space_0);
  // std::cout << "excluding homogeneous_mode_space: \n"
  //           << pretty(dof_space_1.basis) << std::endl;
  EXPECT_EQ(dof_space_1.dim, dof_space_shape_excluding_homogeneous_modes.first);
  EXPECT_EQ(dof_space_1.subspace_dim,
            dof_space_shape_excluding_homogeneous_modes.second);
}

TEST_F(DebugLocalDoFSpaceTest, Test1) {  // failed original method
  // FCC base structure,
  // - with corner atoms {A, B} allowed to displace in 3d
  // - with face atoms {C, D} allowed to displace in the face plane,
  // -> test case where all atoms may displace, but along different dimensions
  //    so no rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");
  Molecule C = Molecule::make_atom("C");
  Molecule D = Molecule::make_atom("D");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),       // AnisoVal type
      {"d0x", "d0y", "d0z"},        // axes names
      Eigen::Matrix3d::Identity(),  // basis
      {}                            // excluded_occs
  };

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1x", "d1y"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d2y", "d2z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d3x", "d3z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.0, lat, FRAC}, {C, D}, {disp_xy}},
       Site{Coordinate{0.0, 0.5, 0.5, lat, FRAC}, {C, D}, {disp_yz}},
       Site{Coordinate{0.5, 0.0, 0.5, lat, FRAC}, {C, D}, {disp_xz}}});

  prim = std::make_shared<config::Prim const>(struc);
  supercell = std::make_shared<config::Supercell const>(
      prim, Eigen::Matrix3l::Identity());

  Index prim_factor_group_size = 48;

  check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
  );

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {9, 9},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {9, 0},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {9, 9},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {9, 9});
}

TEST_F(DebugLocalDoFSpaceTest, Test2) {  // failed original method
  // FCC base structure,
  // - no corner atoms
  // - with face atoms {A, B} allowed to displace in the face plane,
  // -> test case where all atoms may displace, but along different dimensions
  //    so no rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xy{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xy
                       {1., 0., 0.},
                       {0., 1., 0.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_yz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1y", "d1z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in yz
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_xz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d2x", "d2z"},          // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xz
                       {1., 0., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.5, 0.5, 0.0, lat, FRAC}, {A, B}, {disp_xy}},
       Site{Coordinate{0.0, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_yz}},
       Site{Coordinate{0.5, 0.0, 0.5, lat, FRAC}, {A, B}, {disp_xz}}});

  prim = std::make_shared<config::Prim const>(struc);
  supercell = std::make_shared<config::Supercell const>(
      prim, Eigen::Matrix3l::Identity());

  Index prim_factor_group_size = 48;

  check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
  );

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {6, 6},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {6, 0},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {6, 6},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {6, 6});
}

TEST_F(DebugLocalDoFSpaceTest, Test3) {  // passes
  // BCC base structure,
  // - corner atoms {A, B} allowed to displace in xyz
  // - body atoms {A, B} allowed to displace in z,
  // -> test case where all atoms may displace, but along different dimensions
  //    so only 1d (z) rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y", "d0z"},   // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_z{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1z"},                 // axes names
      Eigen::MatrixXd({        // basis: allow displacements in z
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_z}}});

  prim = std::make_shared<config::Prim const>(struc);
  supercell = std::make_shared<config::Supercell const>(
      prim, Eigen::Matrix3l::Identity());

  Index prim_factor_group_size = 16;

  check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
  );

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {4, 4},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {4, 1},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {4, 3},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {4, 3});
}

TEST_F(DebugLocalDoFSpaceTest, Test4) {  // failed for 2x2x2, passed for 2x1x1
  // BCC base structure, non-primitive supercell
  // - corner atoms {A, B} allowed to displace in xyz
  // - body atoms {A, B} allowed to displace in z,
  // -> test case where all atoms may displace, but along different dimensions
  //    so only 1d (z) rigid translations are possible
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d0x", "d0y", "d0z"},   // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  SiteDoFSet disp_z{
      AnisoValTraits::disp(),  // AnisoVal type
      {"d1z"},                 // axes names
      Eigen::MatrixXd({        // basis: allow displacements in z
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}},
       Site{Coordinate{0.5, 0.5, 0.5, lat, FRAC}, {A, B}, {disp_z}}});

  prim = std::make_shared<config::Prim const>(struc);

  Eigen::Matrix3l T;
  T << 2, 0, 0, 0, 2, 0, 0, 0, 2;

  supercell = std::make_shared<config::Supercell const>(prim, T);

  Index prim_factor_group_size = 16;
  Index vol = T.determinant();
  Index prim_disp_dof_space_dim = 4;
  Index homogeneous_mode_dim = 1;

  check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
  );

  check_ExcludeHomogeneousModeSpace(
      // std::pair<Index, Index> initial_dof_space_shape
      {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
      // std::pair<Index, Index> homogeneous_mode_space_shape
      {prim_disp_dof_space_dim * vol, homogeneous_mode_dim},
      // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
      {prim_disp_dof_space_dim * vol,
       prim_disp_dof_space_dim * vol - homogeneous_mode_dim},
      // std::pair<Index, Index> symmetry_adapted_dof_space_shape
      {prim_disp_dof_space_dim * vol,
       prim_disp_dof_space_dim * vol - homogeneous_mode_dim});
}

// previously failed for transformation_matrix_to_super:
//  0  1 -2
//  0  1  2
//  1 -1  0
TEST_F(DebugLocalDoFSpaceTest, Test5) {
  // FCC, 3d displacements, various supercells
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  SiteDoFSet disp_xyz{
      AnisoValTraits::disp(),  // AnisoVal type
      {"dx", "dy", "dz"},      // axes names
      Eigen::MatrixXd({        // basis: allow displacements in xyz
                       {1., 0., 0.},
                       {0., 1., 0.},
                       {0., 0., 1.}})
          .transpose(),
      {}  // excluded_occs
  };

  Lattice lat{Eigen::Vector3d{2.0, 2.0, 0.0}, Eigen::Vector3d{0.0, 2.0, 2.0},
              Eigen::Vector3d{2.0, 0.0, 2.0}};

  BasicStructure struc{lat};
  struc.set_basis(
      {Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {A, B}, {disp_xyz}}});

  prim = std::make_shared<config::Prim const>(struc);

  Index prim_factor_group_size = 48;
  Index prim_disp_dof_space_dim = 3;
  Index homogeneous_mode_dim = 3;

  // lattice enumeration
  int begin_volume = 2;
  int end_volume = 5;
  std::string dirs = "abc";
  Eigen::Matrix3i generating_matrix = Eigen::Matrix3i::Identity();
  ScelEnumProps enumeration_params{begin_volume, end_volume, dirs,
                                   generating_matrix};
  SuperlatticeEnumerator enumerator{prim->basicstructure.lattice(),
                                    prim->sym_info.point_group->element,
                                    enumeration_params};

  // for various supercells:
  for (Lattice const &superlattice : enumerator) {
    Lattice canonical_superlattice = xtal::canonical::equivalent(
        superlattice, prim->sym_info.point_group->element, superlattice.tol());

    supercell =
        std::make_shared<config::Supercell const>(prim, canonical_superlattice);
    Eigen::Matrix3l T =
        supercell->superlattice.transformation_matrix_to_super();

    // std::cout << "--- begin supercell ---" << std::endl;
    // std::cout << "transformation_matrix_to_super:" << std::endl;
    // std::cout << T << std::endl;

    Index vol = T.determinant();

    check_FactorGroupSize(prim_factor_group_size  // Index factor_group_size
    );

    check_ExcludeHomogeneousModeSpace(
        // std::pair<Index, Index> initial_dof_space_shape
        {prim_disp_dof_space_dim * vol, prim_disp_dof_space_dim * vol},
        // std::pair<Index, Index> homogeneous_mode_space_shape
        {prim_disp_dof_space_dim * vol, homogeneous_mode_dim},
        // std::pair<Index, Index> dof_space_shape_excluding_homogeneous_modes
        {prim_disp_dof_space_dim * vol,
         prim_disp_dof_space_dim * vol - homogeneous_mode_dim},
        // std::pair<Index, Index> symmetry_adapted_dof_space_shape
        {prim_disp_dof_space_dim * vol,
         prim_disp_dof_space_dim * vol - homogeneous_mode_dim});
  }
}
