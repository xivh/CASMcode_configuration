#ifndef CASM_config_Configuration
#define CASM_config_Configuration

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace config {

/// \brief Data structure encapsulating configuration DoF values and all
/// information necessary to apply symmetry operations
struct Configuration : public Comparisons<CRTPBase<Configuration>> {
  Configuration(std::shared_ptr<Supercell const> const &_supercell);

  Configuration(std::shared_ptr<Supercell const> const &_supercell,
                clexulator::ConfigDoFValues const &_dof_values);

  std::shared_ptr<Supercell const> supercell;

  clexulator::ConfigDoFValues dof_values;

  /// \brief Less than comparison of Configuration
  bool operator<(Configuration const &rhs) const;

 private:
  friend struct Comparisons<CRTPBase<Configuration>>;

  /// \brief Equality comparison of Configuration
  bool eq_impl(Configuration const &rhs) const;
};

class SupercellSymOp;

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration &apply(SupercellSymOp const &op, Configuration &configuration);

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
/// Configuration
Configuration copy_apply(SupercellSymOp const &op, Configuration configuration);

struct ConfigurationWithProperties {
  explicit ConfigurationWithProperties(
      Configuration const &_configuration,
      std::map<std::string, Eigen::MatrixXd> const &_local_properties = {},
      std::map<std::string, Eigen::VectorXd> const &_global_properties = {});

  Configuration configuration;
  std::map<std::string, Eigen::MatrixXd> local_properties;
  std::map<std::string, Eigen::VectorXd> global_properties;
};

class SupercellSymOp;

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
///     a configuration with properties
ConfigurationWithProperties &apply(SupercellSymOp const &op,
                                   ConfigurationWithProperties &configuration);

/// \brief Apply a symmetry operation specified by a SupercellSymOp to
///     a configuration with properties
ConfigurationWithProperties copy_apply(
    SupercellSymOp const &op, ConfigurationWithProperties configuration);

}  // namespace config
}  // namespace CASM

#endif
