#ifndef CASM_config_DoFSpaceAxisInfo
#define CASM_config_DoFSpaceAxisInfo

#include <optional>
#include <set>
#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}  // namespace xtal

namespace config {

struct DoFSpaceAxisInfo {
  DoFSpaceAxisInfo(
      DoFKey dof_key, xtal::BasicStructure const &prim,
      std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
      std::optional<std::set<Index>> const &sites);

  /// Names the DoF corresponding to each dimension (row) of the DoFSpace basis
  std::vector<std::string> glossary;

  /// \brief The supercell site_index corresponding to each dimension (row) of
  /// the DoFSpace basis. Has value for site DoF only.
  ///
  /// Note:
  /// - For continuous DoF this gives the column index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the site index in `ConfigDoF.occupation()`
  /// corresponding to each row in the DoFSpace basis.
  std::optional<std::vector<Index>> site_index;

  /// \brief The local DoF site DoFSet component index corresponding to each
  /// dimension (row) of the DoFSpace basis. Has value for site DoF only.
  ///
  /// Note:
  /// - For continuous DoF this gives the row index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the index into `Site.occupant_dof()`, which is
  /// the value of `ConfigDoF.occupation()[site_index]`.
  std::optional<std::vector<Index>> dof_component;

  /// Use the site_index and dof_component to lookup the corresponding DoFSpace
  /// basis row index. Has value for site DoF only.
  ///
  /// Usage:
  /// Index basis_row_index = (*m_basis_row_index)[site_index][dof_component]
  std::optional<std::vector<std::vector<Index>>> basis_row_index;
};

}  // namespace config
}  // namespace CASM

#endif
