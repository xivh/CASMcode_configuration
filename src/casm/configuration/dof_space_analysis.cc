#include "casm/configuration/dof_space_analysis.hh"

#include "casm/casm_io/Log.hh"
#include "casm/configuration/DoFSpace_functions.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/group/subgroups.hh"

namespace CASM {
namespace config {

DoFSpaceAnalysisResults::DoFSpaceAnalysisResults(
    clexulator::DoFSpace _symmetry_adapted_dof_space,
    irreps::VectorSpaceSymReport _symmetry_report)
    : symmetry_adapted_dof_space(std::move(_symmetry_adapted_dof_space)),
      symmetry_report(std::move(_symmetry_report)){};

/// \param dof_space_in The DoFSpace for which a symmetry adapted basis is
/// constructed. \param prim The prim \param configuration If null, use the full
/// symmetry of the DoFSpace. If has_value,
///     use the symmetry of the configuration.
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
/// \param calc_wedges If true, calculate the irreducible wedges for the vector
///     space. This may take a long time.
/// \param log Optional logger. If has value and `log->verbosity() >=
/// Log::verbose`,
///     prints step-by-step results to log.
DoFSpaceAnalysisResults dof_space_analysis(
    clexulator::DoFSpace const &dof_space_in, std::shared_ptr<Prim const> prim,
    std::optional<Configuration> configuration,
    std::optional<bool> exclude_homogeneous_modes,
    bool include_default_occ_modes,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ,
    std::optional<std::map<Index, int>> site_index_to_default_occ,
    bool calc_wedges, std::optional<Log> log) {
  if (dof_space_in.basis.cols() == 0) {
    std::stringstream msg;
    msg << "Error in dof_space_analysis: "
        << "Initial DoF space: basis.cols() == 0";
    throw dof_space_analysis_error(msg.str());
  }

  std::shared_ptr<Supercell const> supercell;
  if (configuration.has_value()) {
    supercell = configuration->supercell;
  } else if (dof_space_in.transformation_matrix_to_super.has_value()) {
    supercell = std::make_shared<Supercell const>(
        prim, *dof_space_in.transformation_matrix_to_super);
  } else {
    supercell =
        std::make_shared<Supercell const>(prim, Eigen::Matrix3l::Identity());
  }

  // --- Construct the standard DoF space ---
  clexulator::DoFSpace dof_space_pre1 =
      exclude_homogeneous_mode_space(dof_space_in, exclude_homogeneous_modes);
  if (dof_space_pre1.basis.cols() == 0) {
    std::stringstream msg;
    msg << "Error in dof_space_analysis: "
        << "After excluding homogeneous mode space: basis.cols() == 0";
    throw dof_space_analysis_error(msg.str());
  }

  clexulator::DoFSpace dof_space = exclude_default_occ_modes(
      dof_space_pre1, include_default_occ_modes,
      sublattice_index_to_default_occ, site_index_to_default_occ);
  if (dof_space.basis.cols() == 0) {
    std::stringstream msg;
    msg << "Error in dof_space_analysis: "
        << "After excluding default occ modes: basis.cols() == 0";
    throw dof_space_analysis_error(msg.str());
  }

  // construct symmetry group based on invariance of dof_space and configuration
  std::vector<SupercellSymOp> group(SupercellSymOp::begin(supercell),
                                    SupercellSymOp::end(supercell));

  if (configuration.has_value()) {
    group = make_invariant_subgroup(*configuration, group.begin(), group.end());
    if (group.size() == 0) {
      throw std::runtime_error(
          "Error in dof_space_analysis: config factor group has size==0.");
    }
  }
  if (dof_space.sites.has_value()) {
    group =
        make_invariant_subgroup(*dof_space.sites, group.begin(), group.end());
    if (group.size() == 0) {
      throw std::runtime_error(
          "Error in dof_space_analysis: due to DoFSpace sites, group has "
          "size==0.");
    }
  }

  // get matrix rep and associated SymGroup
  // (for global DoF, this makes the point group, removing duplicates)
  std::shared_ptr<SymGroup const> symgroup;
  std::vector<Eigen::MatrixXd> matrix_rep =
      make_matrix_rep(group, dof_space.dof_key, dof_space.sites, symgroup);

  // use the entire group for irrep decomposition
  std::set<Index> group_indices;
  for (Index i = 0; i < matrix_rep.size(); ++i) {
    group_indices.insert(i);
  }

  // functions to construct sub groups, used to find high symmetry directions
  std::function<irreps::GroupIndicesOrbitSet()> make_cyclic_subgroups_f =
      [=]() { return group::make_cyclic_subgroups(*symgroup); };
  std::function<irreps::GroupIndicesOrbitSet()> make_all_subgroups_f = [=]() {
    return group::make_all_subgroups(*symgroup);
  };

  bool allow_complex = true;

  irreps::IrrepDecomposition irrep_decomposition(
      matrix_rep, group_indices, dof_space.basis, make_cyclic_subgroups_f,
      make_all_subgroups_f, allow_complex, log);

  // Generate report, based on constructed inputs
  irreps::VectorSpaceSymReport symmetry_report = vector_space_sym_report(
      irrep_decomposition, calc_wedges, dof_space.axis_info.glossary);

  // check for error occuring for "disp"
  if (symmetry_report.symmetry_adapted_subspace.cols() <
      dof_space.basis.cols()) {
    std::stringstream msg;
    msg << "Error in dof_space_analysis: "
        << "symmetry_adapted_subspace.cols() < dof_space.basis.cols()";
    throw dof_space_analysis_error(msg.str());
  }

  clexulator::DoFSpace symmetry_adapted_dof_space = clexulator::make_dof_space(
      dof_space.dof_key, dof_space.prim,
      supercell->superlattice.transformation_matrix_to_super(), dof_space.sites,
      symmetry_report.symmetry_adapted_subspace);

  return DoFSpaceAnalysisResults(std::move(symmetry_adapted_dof_space),
                                 std::move(symmetry_report));
}

}  // namespace config
}  // namespace CASM
