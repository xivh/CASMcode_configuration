#include "casm/configuration/Configuration.hh"

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/Prim.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

namespace test {
Eigen::MatrixXd coordinate_frac(
    std::shared_ptr<config::Supercell const> const &supercell) {
  auto const &converter = supercell->unitcellcoord_index_converter;
  auto const &xtal_prim = supercell->prim->basicstructure;
  Index n_sites = converter.total_sites();
  Eigen::MatrixXd R(3, n_sites);
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    R.col(l) = bijk.coordinate(*xtal_prim).const_frac();
  }
  return R;
}

Eigen::MatrixXd coordinate_cart(
    std::shared_ptr<config::Supercell const> const &supercell) {
  auto const &converter = supercell->unitcellcoord_index_converter;
  auto const &xtal_prim = supercell->prim->basicstructure;
  Index n_sites = converter.total_sites();
  Eigen::MatrixXd R(3, n_sites);
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    R.col(l) = bijk.coordinate(*xtal_prim).const_cart();
  }
  return R;
}

}  // namespace test

TEST(ConfigurationTest, Test1) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  std::shared_ptr<config::Supercell const> supercell =
      std::make_shared<config::Supercell const>(prim, T);

  config::Configuration configuration(supercell);

  auto &dof_values = configuration.dof_values;
  EXPECT_EQ(dof_values.occupation.size(), 4);
  EXPECT_EQ(dof_values.local_dof_values.size(), 0);
  EXPECT_EQ(dof_values.global_dof_values.size(), 0);

  Eigen::VectorXi expected(4);
  expected << 0, 0, 0, 0;
  EXPECT_EQ(dof_values.occupation, expected);
}

TEST(ConfigurationJsonTest, Test1) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  std::shared_ptr<config::Supercell const> supercell =
      std::make_shared<config::Supercell const>(prim, T);

  config::Configuration default_configuration(supercell);
  config::Configuration configuration_out = default_configuration;

  auto &dof_values_out = configuration_out.dof_values;
  dof_values_out.occupation(0) = 1;

  jsonParser json;
  to_json(configuration_out, json);
  //  std::cout << json << std::endl;

  {
    config::Configuration configuration_in = default_configuration;
    configuration_in =
        jsonConstructor<config::Configuration>::from_json(json, prim);
    auto &dof_values_in = configuration_in.dof_values;

    EXPECT_EQ(dof_values_in.occupation.size(), 4);
    EXPECT_EQ(dof_values_in.local_dof_values.size(), 0);
    EXPECT_EQ(dof_values_in.global_dof_values.size(), 0);

    Eigen::VectorXi expected(4);
    expected << 1, 0, 0, 0;
    EXPECT_EQ(dof_values_in.occupation, expected);
  }

  {
    std::unique_ptr<config::Configuration> configuration_in;
    configuration_in =
        jsonMake<config::Configuration>::make_from_json(json, prim);
    EXPECT_TRUE(configuration_in != nullptr);

    auto &dof_values_in = configuration_in->dof_values;

    EXPECT_EQ(dof_values_in.occupation.size(), 4);
    EXPECT_EQ(dof_values_in.local_dof_values.size(), 0);
    EXPECT_EQ(dof_values_in.global_dof_values.size(), 0);

    Eigen::VectorXi expected(4);
    expected << 1, 0, 0, 0;
    EXPECT_EQ(dof_values_in.occupation, expected);
  }
}

TEST(ConfigurationWithPropertiesJsonTest, Test1) {
  std::shared_ptr<config::Prim const> prim =
      config::make_shared_prim(test::FCC_binary_prim());

  Eigen::Matrix3l T;
  T << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  std::shared_ptr<config::Supercell const> supercell =
      std::make_shared<config::Supercell const>(prim, T);

  config::ConfigurationWithProperties default_configuration(supercell);
  config::ConfigurationWithProperties configuration_out = default_configuration;

  //  std::cout << "coordinate_cart: \n"
  //            << test::coordinate_cart(supercell).transpose() << std::endl;
  auto &dof_values_out = configuration_out.configuration.dof_values;
  dof_values_out.occupation << 1, 1, 0, 0;

  auto begin = config::SupercellSymOp::begin(supercell);
  auto end = config::SupercellSymOp::end(supercell);
  auto invariant_subgroup =
      make_invariant_subgroup(configuration_out.configuration, begin, end);

  auto &disp = configuration_out.local_properties["disp"];
  disp.resize(3, 4);
  disp.col(0) << 0.01, 0.05, 0.09;
  disp.col(1) << 0.02, 0.06, 0.1;
  disp.col(2) << 0.03, 0.07, 0.11;
  disp.col(3) << 0.04, 0.08, 0.12;
  auto &Hstrain = configuration_out.global_properties["Hstrain"];
  Hstrain.resize(6);
  Hstrain << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06;
  auto &energy = configuration_out.global_properties["energy"];
  energy.resize(1);
  energy << 0.1;

  jsonParser json;
  to_json(configuration_out, json);
  //  std::cout << "~~~ configuration_out ~~~" << std::endl;
  //  std::cout << json << std::endl;

  int i = 0;
  for (auto const &supercell_symop : invariant_subgroup) {
    auto tmp = copy_apply(supercell_symop, configuration_out);
    jsonParser json;
    to_json(tmp, json);
    //    std::cout << "~~~ i=" << i << " ~~~" << std::endl;
    //    std::cout << json << std::endl;
    ++i;
  }

  {
    config::ConfigurationWithProperties configuration_in =
        default_configuration;
    configuration_in =
        jsonConstructor<config::ConfigurationWithProperties>::from_json(json,
                                                                        prim);
    auto &dof_values_in = configuration_in.configuration.dof_values;

    EXPECT_EQ(dof_values_in.occupation.size(), 4);
    EXPECT_EQ(dof_values_in.local_dof_values.size(), 0);
    EXPECT_EQ(dof_values_in.global_dof_values.size(), 0);

    Eigen::VectorXi expected(4);
    expected << 1, 1, 0, 0;
    EXPECT_EQ(dof_values_in.occupation, expected);
  }

  {
    std::unique_ptr<config::ConfigurationWithProperties> configuration_in;
    configuration_in =
        jsonMake<config::ConfigurationWithProperties>::make_from_json(json,
                                                                      prim);
    EXPECT_TRUE(configuration_in != nullptr);

    auto &dof_values_in = configuration_in->configuration.dof_values;

    EXPECT_EQ(dof_values_in.occupation.size(), 4);
    EXPECT_EQ(dof_values_in.local_dof_values.size(), 0);
    EXPECT_EQ(dof_values_in.global_dof_values.size(), 0);

    Eigen::VectorXi expected(4);
    expected << 1, 1, 0, 0;
    EXPECT_EQ(dof_values_in.occupation, expected);
  }
}
