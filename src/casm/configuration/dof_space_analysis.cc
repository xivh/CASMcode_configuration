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
    irreps::IrrepDecomposition _irrep_decomposition,
    irreps::VectorSpaceSymReport _symmetry_report)
    : symmetry_adapted_dof_space(std::move(_symmetry_adapted_dof_space)),
      irrep_decomposition(std::move(_irrep_decomposition)),
      symmetry_report(std::move(_symmetry_report)) {};

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
/// \param symmetrization Method to use for symmetrization of irreducible
///     subspaces. Options are:
///     - "none": Leave the irreducible subspace bases as initially found,
///       reducing computation time.
///     - "fast": Symmetrize the irreducible subspace bases to align along
///       high-symmetry directions using cyclic subgroups. This may not be a
///       complete symmetrization, but is generally fast.
///     - "complete": Symmetrize the irreducible subspace bases to align
///       along high-symmetry directions using all subgroups. For large
///       spaces, finding all subgroups is slow.
/// \param max_iter Maximum number of iterations to use when finding
///     irreducible subspaces. If a non-irreducible subspace cannot be
///     decomposed within this number of iterations, `complete_decomposition`
///     will be set to False. Starting with a different `init_subspace` may
///     result in a complete decomposition.
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
    std::string symmetrization, Index max_iter, bool calc_wedges,
    std::optional<Log> log) {
  if (log.has_value()) {
    log->begin<Log::standard>("DoF space analysis");
    log->indent() << std::endl;
  }
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

  if (log.has_value()) {
    log->indent() << "Initial DoF space dim: " << dof_space_in.basis.cols()
                  << std::endl;
  }

  clexulator::DoFSpace dof_space_pre1 =
      exclude_homogeneous_mode_space(dof_space_in, exclude_homogeneous_modes);
  if (dof_space_pre1.basis.cols() == 0) {
    std::stringstream msg;
    msg << "Error in dof_space_analysis: "
        << "After excluding homogeneous mode space: basis.cols() == 0";
    throw dof_space_analysis_error(msg.str());
  }
  if (log.has_value()) {
    if (dof_space_pre1.basis.cols() != dof_space_in.basis.cols()) {
      log->indent() << "Exclude homogeneous modes." << std::endl;
    }
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
  if (log.has_value()) {
    if (dof_space.basis.cols() != dof_space_pre1.basis.cols()) {
      log->indent() << "Exclude default occupation modes." << std::endl;
    }
    log->indent() << "Final DoF space dim: " << dof_space.basis.cols()
                  << std::endl
                  << std::endl;
  }

  // construct symmetry group based on invariance of dof_space and configuration
  if (log.has_value()) {
    log->custom<Log::standard>("Construct symmetry group");
    log->indent() << std::endl;
    log->indent() << "Add supercell symmetry operations..." << std::endl;
  }

  std::vector<SupercellSymOp> group(SupercellSymOp::begin(supercell),
                                    SupercellSymOp::end(supercell));

  if (configuration.has_value()) {
    if (log.has_value()) {
      log->indent() << "Make configuration invariant subgroup..." << std::endl;
    }
    group = make_invariant_subgroup(*configuration, group.begin(), group.end());
    if (group.size() == 0) {
      throw std::runtime_error(
          "Error in dof_space_analysis: config factor group has size==0.");
    }
  }
  if (dof_space.sites.has_value()) {
    if (log.has_value()) {
      log->indent() << "Make sites invariant subgroup..." << std::endl;
    }
    group =
        make_invariant_subgroup(*dof_space.sites, group.begin(), group.end());
    if (group.size() == 0) {
      throw std::runtime_error(
          "Error in dof_space_analysis: due to DoFSpace sites, group has "
          "size==0.");
    }
  }
  if (log.has_value()) {
    log->indent() << "Number of group elements = " << group.size() << std::endl
                  << std::endl;
  }

  // get matrix rep and associated SymGroup
  // (for global DoF, this makes the point group, removing duplicates)

  if (log.has_value()) {
    log->custom<Log::standard>("Construct matrix representation");
    log->indent() << std::endl;
    if (symmetrization == "none") {
      log->indent() << "Make matrix representation..." << std::endl;
    } else {
      log->indent() << "Make matrix representation, multiplication table, and "
                       "inverse table..."
                    << std::endl;
    }
  }

  std::shared_ptr<SymGroup const> symgroup;
  bool make_symgroup = true;
  if (symmetrization == "none") {
    make_symgroup = false;
  }

  std::vector<Eigen::MatrixXd> matrix_rep = make_matrix_rep(
      group, dof_space.dof_key, dof_space.sites, symgroup, make_symgroup);

  if (symmetrization != "none") {
    log->indent() << "Make Minimal generators... " << std::endl;
    log->indent() << "- Number of elements: " << symgroup->element.size()
                  << std::endl;
    std::set<Index> indices;
    for (Index i = 0; i < symgroup->element.size(); ++i) {
      indices.insert(i);
    }
    group::MakeMinimalSubsetGenerators<xtal::SymOp> x(*symgroup, indices);
    std::set<Index> minimal_generators = x.generators;
    log->indent() << "- Number of generators: " << minimal_generators.size()
                  << std::endl;
    log->indent() << "- Generators: ";
    for (Index j : minimal_generators) {
      log->ostream() << j << " ";
    }
    log->ostream() << std::endl;
  }

  // use the entire group for irrep decomposition
  std::set<Index> group_indices;
  for (Index i = 0; i < matrix_rep.size(); ++i) {
    group_indices.insert(i);
  }

  if (log.has_value()) {
    log->indent() << "Matrix representation: DONE" << std::endl << std::endl;
  }

  // functions to construct sub groups, used to find high symmetry directions
  std::function<irreps::GroupIndicesOrbitSet()> make_cyclic_subgroups_f =
      group::MakeCyclicSubgroups<xtal::SymOp>(symgroup);
  std::function<irreps::GroupIndicesOrbitSet()> make_all_subgroups_f =
      group::MakeAllSubgroups<xtal::SymOp>(symgroup);

  bool allow_complex = true;

  // Note: this is logged internally
  irreps::IrrepDecomposition irrep_decomposition(
      matrix_rep, group_indices, dof_space.basis, make_cyclic_subgroups_f,
      make_all_subgroups_f, allow_complex, symmetrization, max_iter, log);

  // Generate report, based on constructed inputs
  if (log.has_value()) {
    log->custom<Log::standard>("Construct symmetry report");
    log->indent() << std::endl;
  }
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

  if (log.has_value()) {
    log->custom<Log::standard>("Construct symmetry adapted subspace");
    log->indent() << std::endl;
  }
  clexulator::DoFSpace symmetry_adapted_dof_space = clexulator::make_dof_space(
      dof_space.dof_key, dof_space.prim,
      supercell->superlattice.transformation_matrix_to_super(), dof_space.sites,
      symmetry_report.symmetry_adapted_subspace);

  if (log.has_value()) {
    log->end<Log::verbose>("DoF space analysis");
    log->indent() << std::endl;
  }

  return DoFSpaceAnalysisResults(std::move(symmetry_adapted_dof_space),
                                 std::move(irrep_decomposition),
                                 std::move(symmetry_report));
}

}  // namespace config
}  // namespace CASM
