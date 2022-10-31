#include "casm/configuration/Configuration.hh"

#include "casm/clexulator/ConfigDoFValuesTools.hh"
#include "casm/configuration/ConfigCompare.hh"
#include "casm/configuration/ConfigIsEquivalent.hh"
#include "casm/configuration/SupercellSymOp.hh"

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
  if (supercell != rhs.supercell) {
    return false;
  }
  double xtal_tol = supercell->prim->basicstructure->lattice().tol();
  ConfigIsEquivalent equal_to(*this, xtal_tol);
  return equal_to(rhs);
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
