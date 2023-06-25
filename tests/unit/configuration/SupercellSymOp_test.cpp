#include "casm/configuration/SupercellSymOp.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

class SupercellSymOpFCCTest : public testing::Test {
 protected:
  SupercellSymOpFCCTest() {
    std::shared_ptr<config::Prim const> prim =
        config::make_shared_prim(test::FCC_binary_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
  std::shared_ptr<config::SymGroup const> symgroup;
};

TEST_F(SupercellSymOpFCCTest, Test1) {
  EXPECT_EQ(supercell->sym_info.factor_group->element.size(), 48);
  EXPECT_EQ(supercell->sym_info.translation_permutations->size(), 4);
  EXPECT_EQ(supercell->sym_info.factor_group_permutations.size(), 48);
}

TEST_F(SupercellSymOpFCCTest, Test2) {
  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  Index count = 0;
  for (auto it = begin; it != end; ++it) {
    ++count;
  }
  EXPECT_EQ(count, 4 * 48);
}

TEST_F(SupercellSymOpFCCTest, Test3) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
  }
  EXPECT_TRUE(almost_equal(occ_count, Eigen::VectorXi::Constant(size, 48)));
}

TEST_F(SupercellSymOpFCCTest, Test4) {
  config::Configuration configuration(supercell);
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  std::set<Index> site_indices;
  for (Index l = 0; l < n_sites; ++l) {
    site_indices.emplace(l);
  }

  std::vector<Eigen::MatrixXd> occ_matrix_rep = make_local_dof_matrix_rep(
      invariant_subgroup, "occ", site_indices, symgroup);

  // for (Index i=0; i<occ_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << occ_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(occ_matrix_rep.size(), 4 * 48);
}

TEST_F(SupercellSymOpFCCTest, Test5) {
  config::Configuration configuration(supercell);
  configuration.dof_values.occupation(0) = 1;
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  std::set<Index> site_indices;
  for (Index l = 0; l < n_sites; ++l) {
    site_indices.emplace(l);
  }

  std::vector<Eigen::MatrixXd> occ_matrix_rep = make_local_dof_matrix_rep(
      invariant_subgroup, "occ", site_indices, symgroup);

  // for (Index i=0; i<occ_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << occ_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(occ_matrix_rep.size(), 48);
}

class SupercellSymOpFCCTernaryGLStrainDispTest : public testing::Test {
 protected:
  SupercellSymOpFCCTernaryGLStrainDispTest() {
    auto prim =
        config::make_shared_prim(test::FCC_ternary_GLstrain_disp_prim());
    Eigen::Matrix3l T;
    T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
  std::shared_ptr<config::SymGroup const> symgroup;
};

TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  dof_values.local_dof_values.at("disp")(0, 0) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(0) = 0.01;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);
  Eigen::MatrixXd disp_count = Eigen::MatrixXd::Zero(3, size);
  Eigen::VectorXd GLstrain_count = Eigen::VectorXd::Zero(6);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
    disp_count += transformed_dof_values.local_dof_values.at("disp");
    GLstrain_count += transformed_dof_values.global_dof_values.at("GLstrain");
  }
  EXPECT_TRUE(almost_equal(occ_count, Eigen::VectorXi::Constant(size, 48)));
  EXPECT_TRUE(almost_equal(disp_count, Eigen::MatrixXd::Zero(3, size)));

  // positive stretch gets transformed to positive x,y,z
  Eigen::VectorXd expected(6);
  expected << 0.01, 0.01, 0.01, 0.0, 0.0, 0.0;
  expected *= 48.0 * 4.0 / 3.0;
  EXPECT_TRUE(almost_equal(GLstrain_count, expected));

  // shear gets transformed to +/- xy,yz,xz, averages to 0.0
  dof_values.global_dof_values.at("GLstrain") = Eigen::VectorXd::Zero(6);
  dof_values.global_dof_values.at("GLstrain")(4) = 0.01;
  GLstrain_count = Eigen::VectorXd::Zero(6);
  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    GLstrain_count += transformed_dof_values.global_dof_values.at("GLstrain");
  }
  EXPECT_TRUE(almost_equal(GLstrain_count, Eigen::VectorXd::Zero(6)));
}

TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestProduct) {
  // test multiplication table
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin_a = config::SupercellSymOp::begin(supercell);
  auto begin_b = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  dof_values.local_dof_values.at("disp")(0, 0) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(0) = 0.01;

  for (auto it_a = begin_a; it_a != end; ++it_a) {
    for (auto it_b = begin_b; it_b != end; ++it_b) {
      clexulator::ConfigDoFValues ab =
          copy_apply(*it_a, copy_apply(*it_b, dof_values));
      auto product_it = it_a * it_b;
      clexulator::ConfigDoFValues prod = copy_apply(*product_it, dof_values);
      EXPECT_TRUE(almost_equal(ab.occupation, prod.occupation));
      EXPECT_TRUE(almost_equal(ab.local_dof_values.at("disp"),
                               prod.local_dof_values.at("disp")));
      EXPECT_TRUE(almost_equal(ab.global_dof_values.at("GLstrain"),
                               prod.global_dof_values.at("GLstrain")));
    }
  }
}

TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestInverse) {
  // test inverse
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  dof_values.local_dof_values.at("disp")(0, 0) = 1.0;
  dof_values.global_dof_values.at("GLstrain")(0) = 0.01;

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed = copy_apply(*it, dof_values);
    clexulator::ConfigDoFValues check = copy_apply(it.inverse(), transformed);
    EXPECT_TRUE(almost_equal(check.occupation, dof_values.occupation));
    EXPECT_TRUE(almost_equal(check.local_dof_values.at("disp"),
                             dof_values.local_dof_values.at("disp")));
    EXPECT_TRUE(almost_equal(check.global_dof_values.at("GLstrain"),
                             dof_values.global_dof_values.at("GLstrain")));
  }
}

// Test make global matrix rep
TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestGlobalMatrixRep1) {
  config::Configuration configuration(supercell);
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  EXPECT_EQ(invariant_subgroup.size(), 4 * 48);

  std::vector<Eigen::MatrixXd> GLstrain_matrix_rep =
      make_global_dof_matrix_rep(invariant_subgroup, "GLstrain", symgroup);

  // for (Index i=0; i<GLstrain_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << GLstrain_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(GLstrain_matrix_rep.size(), 48);
}

// Test make global matrix rep
TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestGlobalMatrixRep2) {
  config::Configuration configuration(supercell);
  configuration.dof_values.occupation(0) = 1;
  configuration.dof_values.occupation(1) = 1;
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  EXPECT_EQ(invariant_subgroup.size(), 2 * 16);

  std::vector<Eigen::MatrixXd> GLstrain_matrix_rep =
      make_global_dof_matrix_rep(invariant_subgroup, "GLstrain", symgroup);

  // for (Index i=0; i<GLstrain_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << GLstrain_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(GLstrain_matrix_rep.size(), 16);
}

// Test make local matrix rep
TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestLocalMatrixRep1) {
  config::Configuration configuration(supercell);
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  std::set<Index> site_indices;
  for (Index l = 0; l < n_sites; ++l) {
    site_indices.emplace(l);
  }
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  EXPECT_EQ(invariant_subgroup.size(), 4 * 48);

  std::vector<Eigen::MatrixXd> disp_matrix_rep = make_local_dof_matrix_rep(
      invariant_subgroup, "disp", site_indices, symgroup);

  // for (Index i=0; i<disp_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << disp_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(disp_matrix_rep.size(), 4 * 48);
}

// Test make local matrix rep
TEST_F(SupercellSymOpFCCTernaryGLStrainDispTest, TestLocalMatrixRep2) {
  config::Configuration configuration(supercell);
  configuration.dof_values.occupation(0) = 1;
  configuration.dof_values.occupation(1) = 1;
  Index n_sites = supercell->unitcellcoord_index_converter.total_sites();

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  std::set<Index> site_indices;
  for (Index l = 0; l < n_sites; ++l) {
    site_indices.emplace(l);
  }
  auto invariant_subgroup = make_invariant_subgroup(configuration, begin, end);
  EXPECT_EQ(invariant_subgroup.size(), 2 * 16);

  std::vector<Eigen::MatrixXd> disp_matrix_rep = make_local_dof_matrix_rep(
      invariant_subgroup, "disp", site_indices, symgroup);

  // for (Index i=0; i<disp_matrix_rep.size(); ++i) {
  //   std::cout << "i: " << i << "\n" << disp_matrix_rep[i] << std::endl <<
  //   std::endl;
  // }
  EXPECT_EQ(disp_matrix_rep.size(), 2 * 16);
}

class SupercellSymOpSimpleCubicIsingTest : public testing::Test {
 protected:
  SupercellSymOpSimpleCubicIsingTest() {
    auto prim = config::make_shared_prim(test::SimpleCubic_ising_prim());
    Eigen::Matrix3l T = Eigen::Matrix3l::Identity() * 2;
    supercell = std::make_shared<config::Supercell const>(prim, T);
  }

  std::shared_ptr<config::Supercell const> supercell;
  std::shared_ptr<config::SymGroup const> symgroup;
};

TEST_F(SupercellSymOpSimpleCubicIsingTest, Test1) {
  config::Configuration configuration(supercell);
  clexulator::ConfigDoFValues dof_values = configuration.dof_values;
  Index size = dof_values.occupation.size();
  EXPECT_TRUE(almost_equal(dof_values.occupation, Eigen::VectorXi::Zero(size)));

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);

  dof_values.occupation(0) = 1;
  Eigen::VectorXi occ_count = Eigen::VectorXi::Zero(size);

  for (auto it = begin; it != end; ++it) {
    clexulator::ConfigDoFValues transformed_dof_values =
        copy_apply(*it, dof_values);
    occ_count += transformed_dof_values.occupation;
  }
  // Expect 48*8 ops -> 7x down, 1x up; 48*8 ops -> 7x up, 1x down
  EXPECT_TRUE(almost_equal(occ_count,
                           Eigen::VectorXi::Constant(size, 1 * 48 + 7 * 48)));
}
