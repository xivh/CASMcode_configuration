#include "casm/configuration/Configuration.hh"

#include "casm/clexulator/ConfigDoFValuesTools.hh"
#include "casm/clexulator/DoFSpace.hh"
#include "casm/configuration/ConfigCompare.hh"
#include "casm/configuration/ConfigIsEquivalent.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/copy_configuration.hh"

namespace CASM {
namespace config {

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell)
    : supercell(_supercell),
      dof_values(clexulator::make_default_config_dof_values(
          _supercell->prim->basicstructure->basis().size(),
          _supercell->superlattice.size(), _supercell->prim->global_dof_info,
          _supercell->prim->local_dof_info)) {}

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell,
                             clexulator::ConfigDoFValues const &_dof_values)
    : supercell(_supercell), dof_values(_dof_values) {}

/// \brief Less than comparison of Configuration
///
/// - Must have the same Prim
/// - Supercell are compared first, then global DoF, then occupation DoF,
///   then local continuous DoF
bool Configuration::operator<(Configuration const &rhs) const {
  ConfigCompare compare(*this);
  return compare(rhs);
}

/// \brief Equality comparison of Configuration
///
/// - Must have the same Prim
/// - Checks that all DoF are the same, within tolerance
bool Configuration::eq_impl(Configuration const &rhs) const {
  if (*supercell != *rhs.supercell) {
    return false;
  }
  double xtal_tol = supercell->prim->basicstructure->lattice().tol();
  ConfigIsEquivalent equal_to(*this, xtal_tol);
  return equal_to(rhs);
}

/// \brief Convert DoF values into the standard basis
clexulator::ConfigDoFValues make_standard_dof_values(
    Configuration const &config) {
  auto &supercell = config.supercell;
  auto &prim = *supercell->prim;
  return clexulator::to_standard_values(
      config.dof_values, prim.basicstructure->basis().size(),
      supercell->unitcell_index_converter.total_sites(), prim.global_dof_info,
      prim.local_dof_info);
}

/// \brief Set DoF values from other, which is in the standard basis
void set_standard_dof_values(Configuration &config,
                             clexulator::ConfigDoFValues const &other) {
  auto &supercell = config.supercell;
  auto &prim = *supercell->prim;
  config.dof_values = clexulator::from_standard_values(
      other, prim.basicstructure->basis().size(),
      supercell->unitcell_index_converter.total_sites(), prim.global_dof_info,
      prim.local_dof_info);
}

/// \brief Set DoF values associated with a DoFSpace coordinate
///
/// Notes:
/// - All DoF values in the DoFSpace are set to match the given coordinate.
/// - DoF values of other types are not changed.
/// - If the DoFSpace is for a local DoF, only sites included in the
///   DoFSpace are set.
/// - This method does not support occupation DoFSpace.
///
/// \param config The configuration being assigned to
/// \param dof_space The DoFSpace. For local DoF spaces, the dof_space must
///     tile into the config supercell.
/// \param dof_space_coordinate The normal coordinate specifying the DoF
///     values in the DoFSpace basis.
void set_dof_space_values(config::Configuration &config,
                          clexulator::DoFSpace const &dof_space,
                          Eigen::VectorXd const &dof_space_coordinate) {
  if (dof_space.dof_key == "occ") {
    std::stringstream msg;
    msg << "Error: CASM::config::set_dof_space_values is not "
           "supported for occupation."
        << std::endl;
    throw std::runtime_error(msg.str());
  }

  // For global DoF, or local DoF if the supercells are identical,
  // can simply directly assign DoF values
  if (dof_space.is_global ||
      (config.supercell->superlattice.transformation_matrix_to_super() ==
       dof_space.transformation_matrix_to_super.value())) {
    clexulator::set_dof_value(
        config.dof_values,
        config.supercell->superlattice.transformation_matrix_to_super(),
        dof_space, dof_space_coordinate);
    return;
  }

  // Otherwise, need to:
  // 1) assign DoF values to the DoF space supercell,
  // 2) make a super configuration in the config supercell,
  // 3) copy DoF space values (will be local DoF only)

  // make dof_space supercell
  auto dof_space_supercell = std::make_shared<Supercell const>(
      config.supercell->prim, dof_space.transformation_matrix_to_super.value());

  // if dof_space supercell does not tile `config.supercell`, throw
  if (config.supercell != dof_space_supercell) {
    auto const &superlattice = config.supercell->superlattice.superlattice();
    auto const &unit_lattice = dof_space_supercell->superlattice.superlattice();
    double tol = std::max(superlattice.tol(), unit_lattice.tol());
    if (!xtal::is_superlattice(superlattice, unit_lattice, tol).first) {
      throw std::runtime_error(
          "Error in Configuration.set_dof_space_values: "
          "configuration is not tiled by the dof_space supercell");
    }
  }

  // assign DoF values to the DoF space supercell
  config::Configuration dof_space_config(dof_space_supercell);
  clexulator::set_dof_value(
      dof_space_config.dof_values,
      dof_space_config.supercell->superlattice.transformation_matrix_to_super(),
      dof_space, dof_space_coordinate);

  // make a super configuration in the config supercell
  config::Configuration dof_space_superconfig =
      copy_configuration(dof_space_config, config.supercell);

  // copy DoF space values (will be local DoF only)
  Eigen::MatrixXd const &x =
      dof_space_superconfig.dof_values.local_dof_values.at(dof_space.dof_key);
  Eigen::MatrixXd &y = config.dof_values.local_dof_values.at(dof_space.dof_key);
  for (Index l : *dof_space.sites) {
    y.col(l) = x.col(l);
  }
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration &apply(SupercellSymOp const &op, Configuration &configuration) {
  apply(op, configuration.dof_values);
  return configuration;
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration copy_apply(SupercellSymOp const &op,
                         Configuration configuration) {
  return apply(op, configuration);
}

ConfigurationWithProperties::ConfigurationWithProperties(
    Configuration const &_configuration,
    std::map<std::string, Eigen::MatrixXd> const &_local_properties,
    std::map<std::string, Eigen::VectorXd> const &_global_properties)
    : configuration(_configuration),
      local_properties(_local_properties),
      global_properties(_global_properties) {}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
///     a configuration with properties
ConfigurationWithProperties &apply(
    SupercellSymOp const &op,
    ConfigurationWithProperties &config_with_properties) {
  apply(op, config_with_properties.configuration);

  Configuration &configuration = config_with_properties.configuration;
  Index n_sites =
      configuration.supercell->unitcellcoord_index_converter.total_sites();
  Index prim_fg_index = op.prim_factor_group_index();
  xtal::SymOp symop = op.to_symop();

  // transform global properties
  for (auto &property : config_with_properties.global_properties) {
    AnisoValTraits traits(property.first);
    Eigen::MatrixXd M = traits.symop_to_matrix(
        get_matrix(symop), get_translation(symop), get_time_reversal(symop));
    property.second = M * property.second;
  }

  // transform then permute local properties
  sym_info::Permutation combined_permute{op.combined_permute()};
  for (auto &property : config_with_properties.local_properties) {
    AnisoValTraits traits(property.first);
    Eigen::MatrixXd M = traits.symop_to_matrix(
        get_matrix(symop), get_translation(symop), get_time_reversal(symop));
    Eigen::MatrixXd tmp = M * property.second;
    // permute values amongst sites
    for (Index l = 0; l < n_sites; ++l) {
      property.second.col(l) = tmp.col(combined_permute[l]);
    }
  }

  return config_with_properties;
}

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
///     a configuration with properties
ConfigurationWithProperties copy_apply(
    SupercellSymOp const &op,
    ConfigurationWithProperties config_with_properties) {
  apply(op, config_with_properties);
  return config_with_properties;
}

}  // namespace config
}  // namespace CASM
