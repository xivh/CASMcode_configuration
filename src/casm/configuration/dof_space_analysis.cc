#include "casm/configuration/dof_space_analysis.hh"

#include "casm/casm_io/Log.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/group/subgroups.hh"

namespace CASM {
namespace config {

// projector(Eigen::MatrixXd::Zero(standard_dof_space.basis.cols(),
//                                 standard_dof_space.basis.cols()))

namespace {  // (anonymous)

}  // namespace

DoFSpaceAnalysisResults::DoFSpaceAnalysisResults(
    clexulator::DoFSpace _symmetry_adapted_dof_space,
    irreps::VectorSpaceSymReport _symmetry_report)
    : symmetry_adapted_dof_space(std::move(_symmetry_adapted_dof_space)),
      symmetry_report(std::move(_symmetry_report)){};

/// \param dof_space Names of degree of freedoms for which the analysis
///     is run. The default includes all DoF types in the prim.
/// \param supercell Map of identifier string -> Configuration
/// \param group Exclude homogeneous modes if this
///     is true, or include if this is false. If this is null (default),
///     exclude homogeneous modes for dof==\"disp\" only.
/// \param include_default_occ_modes Include the dof component for the
///     default occupation value on each site with occupation DoF. The
///     default is to exclude these modes because they are not
///     independent. This parameter is only checked dof==\"occ\".
DoFSpaceAnalysisResults dof_space_analysis(
    clexulator::DoFSpace const &dof_space, std::shared_ptr<Prim const> prim,
    std::optional<Configuration> configuration, bool calc_wedges) {
  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  if (dof_space.transformation_matrix_to_super.has_value()) {
    T = *dof_space.transformation_matrix_to_super;
  }
  auto supercell = std::make_shared<Supercell const>(prim, T);

  // construct symmetry group based on invariance of dof_space and configuration
  std::vector<SupercellSymOp> group(SupercellSymOp::begin(supercell),
                                    SupercellSymOp::end(supercell));

  if (configuration.has_value()) {
    group = make_invariant_subgroup(*configuration, group.begin(), group.end());
  }
  if (dof_space.sites.has_value()) {
    group =
        make_invariant_subgroup(*dof_space.sites, group.begin(), group.end());
  }

  // get matrix rep and associated SymGroup
  // (for global DoF, this makes the point group, removing duplicates)
  std::shared_ptr<SymGroup const> symgroup;
  std::vector<Eigen::MatrixXd> matrix_rep =
      make_matrix_rep(group, dof_space.dof_key, dof_space.sites, symgroup);

  // use the entire group for irrep decomposition
  std::set<Index> group_indices;
  for (Index i = 0; i < group.size(); ++i) {
    group_indices.insert(i);
  }

  // functions to construct sub groups, used to find high symmetry directions
  std::function<irreps::GroupIndicesOrbitSet()> make_cyclic_subgroups_f =
      [=]() { return group::make_cyclic_subgroups(*symgroup); };
  std::function<irreps::GroupIndicesOrbitSet()> make_all_subgroups_f = [=]() {
    return group::make_all_subgroups(*symgroup);
  };

  bool allow_complex = true;

  irreps::VectorSpaceSymReport symmetry_report;
  try {
    irreps::IrrepDecomposition irrep_decomposition(
        matrix_rep, group_indices, dof_space.basis, make_cyclic_subgroups_f,
        make_all_subgroups_f, allow_complex);

    // Generate report, based on constructed inputs
    symmetry_report = vector_space_sym_report(irrep_decomposition, calc_wedges);

  } catch (std::exception &e) {
    CASM::err_log() << "Error in dof_space_analysis: " << e.what() << std::endl;
    throw e;
  }
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
