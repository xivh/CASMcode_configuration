#include "casm/configuration/DoFSpace_functions.hh"

namespace CASM {
namespace config {

/// Removes specified occupation modes from the DoFSpace basis, by sublattice
///
/// \param dof_space Initial DoF space
/// \param sublattice_index_to_default_occ Table of prim sublattice index to
/// occupation
///     index to treat as the default occupant and remove
/// \return DoFSpace with basis updated by setting rows corresponding to the
/// default
///     occupant to zero, and then columns that are entirely zero are removed.
clexulator::DoFSpace exclude_default_occ_modes_by_sublattice(
    clexulator::DoFSpace const &dof_space,
    std::map<int, int> sublattice_index_to_default_occ) {
  if (dof_space.dof_key != "occ") {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_sublattice: Not occupation DoF");
  }
  if (!dof_space.transformation_matrix_to_super.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_sublattice: no "
        "transformation_matrix_to_super");
  }
  if (!dof_space.axis_info.site_index.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_sublattice: no site_index");
  }
  if (!dof_space.axis_info.dof_component.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_sublattice: no dof_component");
  }

  xtal::BasicStructure const &prim = *dof_space.prim;
  auto const &T = *dof_space.transformation_matrix_to_super;
  xtal::UnitCellCoordIndexConverter l_to_bijk(T, prim.basis().size());

  // validate sublattice_index_to_default_occ
  for (auto const &pair : sublattice_index_to_default_occ) {
    int b = pair.first;
    int default_occ = pair.second;
    if (b < 0 || b >= prim.basis().size()) {
      std::stringstream ss;
      ss << "Error in exclude_default_occ_modes_by_sublattice: sublattice=" << b
         << " is out of range" << std::endl;
      throw std::runtime_error(ss.str());
    }
    auto const &site = prim.basis()[b];
    if (default_occ < 0 || default_occ >= site.occupant_dof().size()) {
      std::stringstream ss;
      ss << "Error in exclude_default_occ_modes_by_sublattice: default_occ="
         << default_occ << " is out of range for sublattice=" << b << std::endl;
      throw std::runtime_error(ss.str());
    }
  }

  // Copy current basis
  Eigen::MatrixXd basis = dof_space.basis;

  // Set rows to zero which correspond to default occ
  std::vector<Index> const &site_index = *dof_space.axis_info.site_index;
  std::vector<Index> const &occ_index = *dof_space.axis_info.dof_component;
  for (Index i = 0; i < basis.rows(); ++i) {
    Index l = site_index[i];
    int b = l_to_bijk(l).sublattice();
    auto it = sublattice_index_to_default_occ.find(b);
    if (it == sublattice_index_to_default_occ.end()) {
      continue;
    }
    int default_occ = it->second;
    if (occ_index[i] == default_occ) {
      for (Index j = 0; j < basis.cols(); ++j) {
        basis(i, j) = 0.0;
      }
    }
  }

  // Copy non-zero columns
  Eigen::MatrixXd tbasis(basis.rows(), basis.cols());
  Index non_zero_cols = 0;
  for (Index j = 0; j < basis.cols(); ++j) {
    if (!almost_zero(basis.col(j))) {
      tbasis.col(non_zero_cols) = basis.col(j);
      ++non_zero_cols;
    }
  }

  // Construct with only non-zero columns
  return clexulator::make_dof_space(dof_space.dof_key, dof_space.prim,
                                    dof_space.transformation_matrix_to_super,
                                    dof_space.sites,
                                    tbasis.leftCols(non_zero_cols));
}

/// Removes specified occupation modes from the DoFSpace basis, by supercell
/// site index
///
/// \param dof_space Initial DoF space
/// \param site_index_to_default_occ Table of supercell site index to occupation
///     index to treat as the default occupant and remove
/// \return DoFSpace with basis updated by setting rows corresponding to the
/// default
///     occupant to zero, and then columns that are entirely zero are removed.
clexulator::DoFSpace exclude_default_occ_modes_by_site(
    clexulator::DoFSpace const &dof_space,
    std::map<Index, int> site_index_to_default_occ) {
  if (dof_space.dof_key != "occ") {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_site: Not occupation DoF");
  }
  if (!dof_space.transformation_matrix_to_super.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_site: no "
        "transformation_matrix_to_super");
  }
  if (!dof_space.axis_info.site_index.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_site: no site_index");
  }
  if (!dof_space.axis_info.dof_component.has_value()) {
    throw std::runtime_error(
        "Error in exclude_default_occ_modes_by_site: no dof_component");
  }

  xtal::BasicStructure const &prim = *dof_space.prim;
  auto const &T = *dof_space.transformation_matrix_to_super;
  xtal::UnitCellCoordIndexConverter l_to_bijk(T, prim.basis().size());

  // validate site_index_to_default_occ
  for (auto const &pair : site_index_to_default_occ) {
    Index l = pair.first;
    if (l < 0 || l >= l_to_bijk.total_sites()) {
      std::stringstream ss;
      ss << "Error in exclude_default_occ_modes_by_site: site=" << l
         << " is out of range" << std::endl;
      throw std::runtime_error(ss.str());
    }
    int b = l_to_bijk(l).sublattice();
    int default_occ = pair.second;
    auto const &site = prim.basis()[b];
    if (default_occ < 0 || default_occ >= site.occupant_dof().size()) {
      std::stringstream ss;
      ss << "Error in exclude_default_occ_modes_by_site: default_occ="
         << default_occ << " is out of range for site=" << l
         << " (sublattice=" << b << ")" << std::endl;
      throw std::runtime_error(ss.str());
    }
  }

  // Copy current basis
  Eigen::MatrixXd basis = dof_space.basis;

  // Set rows to zero which correspond to default occ
  std::vector<Index> const &site_index = *dof_space.axis_info.site_index;
  std::vector<Index> const &occ_index = *dof_space.axis_info.dof_component;
  for (Index i = 0; i < basis.rows(); ++i) {
    Index l = site_index[i];
    auto it = site_index_to_default_occ.find(l);
    if (it == site_index_to_default_occ.end()) {
      continue;
    }
    int default_occ = it->second;
    if (occ_index[i] == default_occ) {
      for (Index j = 0; j < basis.cols(); ++j) {
        basis(i, j) = 0.0;
      }
    }
  }

  // Copy non-zero columns
  Eigen::MatrixXd tbasis(basis.rows(), basis.cols());
  Index non_zero_cols = 0;
  for (Index j = 0; j < basis.cols(); ++j) {
    if (!almost_zero(basis.col(j))) {
      tbasis.col(non_zero_cols) = basis.col(j);
      ++non_zero_cols;
    }
  }

  // Construct with only non-zero columns
  return clexulator::make_dof_space(dof_space.dof_key, dof_space.prim,
                                    dof_space.transformation_matrix_to_super,
                                    dof_space.sites,
                                    tbasis.leftCols(non_zero_cols));
}

/// Removes specified occupation modes from the DoFSpace basis
///
/// \param dof_space_in The DoFSpace for which a symmetry adapted basis is
///     constructed.
/// \param exclude_homogeneous_modes Exclude homogeneous modes if this
///     is true, or include if this is false. If this is null (default),
///     exclude homogeneous modes for dof==\"disp\" only.
/// \param include_default_occ_modes Include the dof component for the
///     default occupation value on each site with occupation DoF. The
///     default is to exclude these modes because they are not
///     independent. This parameter is only checked dof==\"occ\". If
///     false, the default occupation is determined using
///     `site_index_to_default_occ` if that is provided, else using
///     `sublattice_index_to_default_occ` if that is provided, else using
///     occupation index 0.
/// \param sublattice_index_to_default_occ Optional values of default
///     occupation index (value), specified by sublattice index (key).
/// \param site_index_to_default_occ Optional values of default
///     occupation index (value), specified by supercell site index (key).

clexulator::DoFSpace exclude_default_occ_modes(
    clexulator::DoFSpace const &dof_space_in, bool include_default_occ_modes,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ,
    std::optional<std::map<Index, int>> site_index_to_default_occ) {
  if (dof_space_in.dof_key == "occ" && !include_default_occ_modes) {
    if (site_index_to_default_occ.has_value()) {
      return exclude_default_occ_modes_by_site(dof_space_in,
                                               *site_index_to_default_occ);
    } else if (sublattice_index_to_default_occ.has_value()) {
      return exclude_default_occ_modes_by_sublattice(
          dof_space_in, *sublattice_index_to_default_occ);
    } else {
      return clexulator::exclude_default_occ_modes(dof_space_in);
    }
  }
  return dof_space_in;
}

/// Removes homogeneous modes from the DoFSpace basis
///
/// \param dof_space_in The DoFSpace for which a symmetry adapted basis is
///     constructed.
/// \param exclude_homogeneous_modes Exclude homogeneous modes if this
///     is true, or include if this is false. If this is null (default),
///     exclude homogeneous modes for dof==\"disp\" only.
clexulator::DoFSpace exclude_homogeneous_mode_space(
    clexulator::DoFSpace const &dof_space_in,
    std::optional<bool> exclude_homogeneous_modes) {
  if (!exclude_homogeneous_modes.has_value()) {
    if (dof_space_in.dof_key == "disp") {
      exclude_homogeneous_modes = true;
    } else {
      exclude_homogeneous_modes = false;
    }
  }

  if (*exclude_homogeneous_modes) {
    return clexulator::exclude_homogeneous_mode_space(dof_space_in);
  } else {
    return dof_space_in;
  }
}

}  // namespace config
}  // namespace CASM
