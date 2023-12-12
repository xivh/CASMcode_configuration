#ifndef CASM_config_config_space_analysis
#define CASM_config_config_space_analysis

#include "casm/clexulator/DoFSpace.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

struct ConfigSpaceAnalysisResults {
  ConfigSpaceAnalysisResults(
      clexulator::DoFSpace const &_standard_dof_space,
      std::map<std::string, std::vector<Eigen::VectorXd>>
          _equivalent_dof_values,
      std::map<std::string, std::vector<Configuration>>
          _equivalent_configurations,
      Eigen::MatrixXd const &_projector, Eigen::VectorXd const &_eigenvalues,
      clexulator::DoFSpace const &_symmetry_adapted_config_space);

  /// \brief Standard DoF space, may exclude default occupation modes
  ///     or homogeneous displacement modes, depending on method options
  clexulator::DoFSpace const standard_dof_space;

  /// \brief DoF values of all equivalent configurations in the
  ///     fully commensurate supercell, expressed in the basis of the
  ///     standard DoF space, with key == input configuration identifier
  std::map<std::string, std::vector<Eigen::VectorXd>> const
      equivalent_dof_values;

  /// \brief All equivalent configurations in the fully commensurate
  ///     supercell, with key == input configuration identifier
  std::map<std::string, std::vector<Configuration>> const
      equivalent_configurations;

  /// \brief Projection matrix
  Eigen::MatrixXd const projector;

  /// \brief Non-zero eigenvalues of projector
  Eigen::VectorXd const eigenvalues;

  /// \brief Symmetry-adapted config space, with basis formed by
  ///     eigenvectors of P
  clexulator::DoFSpace const symmetry_adapted_dof_space;
};

std::map<DoFKey, ConfigSpaceAnalysisResults> config_space_analysis(
    std::map<std::string, Configuration> const &configurations,
    std::optional<std::vector<DoFKey>> dofs = std::nullopt,
    std::optional<bool> exclude_homogeneous_modes = std::nullopt,
    bool include_default_occ_modes = false,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ =
        std::nullopt,
    std::optional<std::map<Index, int>> site_index_to_default_occ =
        std::nullopt,
    double tol = TOL);

}  // namespace config
}  // namespace CASM

#endif
