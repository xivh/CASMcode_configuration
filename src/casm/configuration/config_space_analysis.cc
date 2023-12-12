#include "casm/configuration/config_space_analysis.hh"

#include "casm/configuration/DoFSpace_functions.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/copy_configuration.hh"
#include "casm/crystallography/CanonicalForm.hh"

namespace CASM {
namespace config {

ConfigSpaceAnalysisResults::ConfigSpaceAnalysisResults(
    clexulator::DoFSpace const &_standard_dof_space,
    std::map<std::string, std::vector<Eigen::VectorXd>> _equivalent_dof_values,
    std::map<std::string, std::vector<Configuration>>
        _equivalent_configurations,
    Eigen::MatrixXd const &_projector, Eigen::VectorXd const &_eigenvalues,
    clexulator::DoFSpace const &_symmetry_adapted_dof_space)
    : standard_dof_space(_standard_dof_space),
      equivalent_dof_values(std::move(_equivalent_dof_values)),
      equivalent_configurations(std::move(_equivalent_configurations)),
      projector(_projector),
      eigenvalues(_eigenvalues),
      symmetry_adapted_dof_space(_symmetry_adapted_dof_space) {}

/// \brief Construct symmetry adapted bases in the DoF space spanned by the
///    set of configurations symmetrically equivalent to the input
///    configurations
///
/// This method:
///
/// 1. Constructs a projector onto the DoF space spanned by all
///    configurations symmetrically equivalent to a set of input
///    configurations, and
/// 2. Finds the eigenvectors spanning that space. The eigenvectors
///    are axes that can be used as a basis for symmetry adapted order
///    parameters.
///
/// This method is faster than the function `dof_space_analysis`, and
/// often the resulting axes do lie along high-symmetry directions, but
/// they may not lie within a single irreducible subspace, and they are
/// not explicitly rotated to the highest symmetry directions.
///
/// \param dofs Names of degree of freedoms for which the analysis
///     is run. The default includes all DoF types in the prim.
/// \param configurations Map of identifier string -> Configuration
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
/// \param tol Tolerance used for identifying zero-valued eigenvalues.
///
/// \returns Results, including project, eigenvalues, and symmetry
///     adapted basis, for each requested DoF type.
std::map<DoFKey, ConfigSpaceAnalysisResults> config_space_analysis(
    std::map<std::string, Configuration> const &configurations,
    std::optional<std::vector<DoFKey>> dofs,
    std::optional<bool> exclude_homogeneous_modes,
    bool include_default_occ_modes,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ,
    std::optional<std::map<Index, int>> site_index_to_default_occ, double tol) {
  std::map<DoFKey, ConfigSpaceAnalysisResults> results;

  if (configurations.size() == 0) {
    return results;
  }
  std::shared_ptr<Prim const> prim =
      configurations.begin()->second.supercell->prim;
  if (!dofs.has_value()) {
    dofs = all_dof_types(*prim->basicstructure);
  }

  // --- Generate primitive configurations ---
  // prim config -> ID (might be duplicates from input configurations)
  std::map<Configuration, std::string> prim_configs;
  for (auto const &pair : configurations) {
    prim_configs.emplace(
        make_in_canonical_supercell(make_primitive(pair.second)), pair.first);
  }

  // --- Generate fully commensurate supercell ---
  std::set<xtal::Lattice> lattices;
  for (auto const &pair : prim_configs) {
    lattices.insert(pair.first.supercell->superlattice.superlattice());
  }
  auto const &fg = prim->sym_info.factor_group->element;
  xtal::Lattice super_lat = xtal::make_fully_commensurate_superduperlattice(
      lattices.begin(), lattices.end(), fg.begin(), fg.end());
  auto const &pg = prim->sym_info.point_group->element;
  super_lat = xtal::canonical::equivalent(super_lat, pg);
  auto shared_supercell = std::make_shared<Supercell const>(prim, super_lat);

  // --- Generate symmetry adapted config spaces ---
  for (auto const &dof_key : *dofs) {
    // --- Construct the standard DoF space ---
    clexulator::DoFSpace dof_space_pre2 = clexulator::make_dof_space(
        dof_key, prim->basicstructure,
        shared_supercell->superlattice.transformation_matrix_to_super());

    clexulator::DoFSpace dof_space_pre1 = exclude_homogeneous_mode_space(
        dof_space_pre2, exclude_homogeneous_modes);

    clexulator::DoFSpace standard_dof_space = exclude_default_occ_modes(
        dof_space_pre1, include_default_occ_modes,
        sublattice_index_to_default_occ, site_index_to_default_occ);

    // --- Begin projector construction ---
    std::map<std::string, std::vector<Eigen::VectorXd>> equivalent_dof_values;
    std::map<std::string, std::vector<Configuration>> equivalent_configurations;
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(standard_dof_space.basis.cols(),
                                              standard_dof_space.basis.cols());

    // std::cout << "P:\n" << P << std::endl;
    for (auto const &prim_config : prim_configs) {
      Configuration prototype =
          copy_configuration(prim_config.first, shared_supercell);
      std::vector<Configuration> equivalents =
          make_equivalents(prototype, SupercellSymOp::begin(shared_supercell),
                           SupercellSymOp::end(shared_supercell));

      std::vector<Eigen::VectorXd> equiv_x;
      for (auto const &config : equivalents) {
        Eigen::VectorXd x = get_normal_coordinate(
            config.dof_values,
            shared_supercell->superlattice.transformation_matrix_to_super(),
            standard_dof_space);

        // clean up x?
        for (int i = 0; i < x.size(); ++i) {
          if (almost_zero(x(i), tol)) {
            x(i) = 0.0;
          }
        }

        P += x * x.transpose();
        equiv_x.push_back(x);
        // std::cout << "P:\n" << P << std::endl;
      }
      equivalent_dof_values[prim_config.second] = equiv_x;
      equivalent_configurations[prim_config.second] = equivalents;
    }
    // std::cout << "P:\n" << P << std::endl;

    // clean up P?
    for (int i = 0; i < P.rows(); ++i) {
      for (int j = 0; j < P.cols(); ++j) {
        if (almost_zero(P(i, j), tol)) {
          P(i, j) = 0.0;
        }
      }
    }

    // --- Eigendecomposition of P ---
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(P);
    Eigen::VectorXd D = solver.eigenvalues();
    Eigen::MatrixXd V = solver.eigenvectors();

    // std::cout << "eigenvalues:\n" << D << std::endl;
    // std::cout << "eigenvectors:\n" << V << std::endl;
    // std::cout << "check:\n" << V * D.asDiagonal() * V.transpose() <<
    // std::endl;

    // --- Identify non-zero eigenvalues and corresponding eigenvectors ---
    Eigen::VectorXd D_nonzero(D.size());
    Eigen::MatrixXd V_nonzero(P.rows(), P.cols());

    int i_nonzero = 0;
    for (int i = 0; i < D.size(); ++i) {
      if (!almost_zero(D(i), tol)) {
        D_nonzero(i_nonzero) = D(i);
        V_nonzero.col(i_nonzero) = V.col(i);
        ++i_nonzero;
      }
    }

    // std::cout << "non-zero eigenvalues:\n"
    //           << D_nonzero.head(i_nonzero) << std::endl;
    // std::cout << "non-zero eigenvectors:\n"
    //           << V_nonzero.leftCols(i_nonzero) << std::endl;
    // std::cout << "B * non-zero eigenvectors:\n"
    //           << standard_dof_space.basis() * V_nonzero.leftCols(i_nonzero)
    //           << std::endl;

    if (i_nonzero == 0) {
      throw std::runtime_error(
          "Error in config_space_analysis: symmetry adapted config space is "
          "null");
    }

    // --- Store results ---
    Eigen::VectorXd eigenvalues = D_nonzero.head(i_nonzero);

    clexulator::DoFSpace symmetry_adapted_dof_space =
        clexulator::make_dof_space(
            standard_dof_space.dof_key, standard_dof_space.prim,
            standard_dof_space.transformation_matrix_to_super,
            standard_dof_space.sites,
            standard_dof_space.basis * V_nonzero.leftCols(i_nonzero));

    results.emplace(std::piecewise_construct, std::forward_as_tuple(dof_key),
                    std::forward_as_tuple(
                        standard_dof_space, std::move(equivalent_dof_values),
                        std::move(equivalent_configurations), P, eigenvalues,
                        symmetry_adapted_dof_space));
  }

  return results;
}

}  // namespace config
}  // namespace CASM
