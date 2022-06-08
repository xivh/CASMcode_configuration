
// Eigen::VectorXd const &delta_prim_dof_values(Index config_site_index, Index
// new_occ, Index old_occ) Eigen::VectorXd const &delta_prim_dof_values(Index
// dof_space_site_index, Index new_occ, Index old_occ) Eigen::VectorXd const
// &delta_order_parameter(Index config_site_index, Index new_occ, Index old_occ)
// Eigen::VectorXd const &delta_order_parameter(Index dof_space_site_index,
// Index new_occ, Index old_occ)

#ifndef CASM_config_DoFSpaceCoordinate
#define CASM_config_DoFSpaceCoordinate

#include <optional>

#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/DoFDecl.hh"

namespace CASM {
namespace config {

/// \brief Method to calculate coordinates & changes in coordinates in a
/// DoFSpace
class DoFSpaceCoordinate {
  /// \brief Constructor - calculate coordinates in a DoFSpace
  DoFSpaceCoordinate(DoFSpace const &dof_space,
                     ConfigDoFValues const *_dof_values = nullptr);

  /// \brief Set internal pointer to DoF values - must be consistent with
  ///     supercell_neighbor_list
  void set(ConfigDoFValues const *_dof_values);

  /// \brief Get internal pointer to DoF values
  ConfigDoFValues const *get() const;

 private:
  /// Holds temporary correlations used in calculations internally
  Eigen::VectorXd m_coordinate;

  /// Holds last delta coordinate results
  Eigen::VectorXd m_delta_coordinate;

  /// Configuration values to use
  Configuration const *m_configuration;

  /// DoFSpace to use
  DoFSpace const m_dof_space;

  /// Pseudo-inverse of m_dof_space.basis
  Eigen::MatrixXd m_basis_pinv;
};

};

}  // namespace config
}  // namespace CASM

#endif
