#include "casm/configuration/irreps/IrrepDecomposition.hh"

#include <random>

#include "casm/configuration/irreps/IrrepDecompositionImpl.hh"
#include "casm/configuration/irreps/VectorSymCompare_v2.hh"
#include "casm/configuration/misc.hh"
#include "casm/misc/CASM_Eigen_math.hh"

// logging
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"
#include "casm/configuration/irreps/misc.hh"

namespace CASM {

namespace irreps {

namespace {
template <int _required_verbosity = Log::standard>
void print_irreps(Log &log, std::string what,
                  std::vector<IrrepInfo> const &irreps) {
  log.begin_section<_required_verbosity>();
  if (log.print()) {
    log.indent() << what << ": " << std::endl;

    jsonParser tmp;
    to_json(irreps, tmp);
    std::stringstream ss;
    ss << tmp << std::endl;
    log.verbatim(ss.str(), true);
    log.indent() << std::endl;
  }
  log.end_section();
}

}  // namespace

IrrepInfo::IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters)
    : trans_mat(std::move(_trans_mat)),
      irrep_dim(trans_mat.rows()),
      vector_dim(trans_mat.cols()),
      characters(std::move(_characters)),
      complex(!almost_zero(trans_mat.imag())),
      pseudo_irrep(false),
      index(0) {}

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXcd const &trans_mat) {
  Eigen::VectorXcd tchar(1);
  tchar(0) = std::complex<double>(double(trans_mat.rows()), 0.);
  return IrrepInfo(trans_mat, tchar);
}

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXd const &trans_mat) {
  Eigen::VectorXcd tchar(1);
  tchar(0) = std::complex<double>(double(trans_mat.rows()), 0.);
  return IrrepInfo(trans_mat.template cast<std::complex<double>>(), tchar);
}

/// \brief Assumes that irreps are real, and concatenates their individual
/// trans_mats to form larger trans_mat
Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps,
                               bool allow_complex) {
  Index row = 0;
  Index col = 0;
  for (auto const &irrep : irreps) {
    col = irrep.vector_dim;
    row += irrep.irrep_dim;
  }

  Eigen::MatrixXd trans_mat(row, col);
  row = 0;

  if (!allow_complex) {
    for (auto const &irrep : irreps) {
      trans_mat.block(row, 0, irrep.irrep_dim, irrep.vector_dim) =
          irrep.trans_mat.real();
      row += irrep.irrep_dim;
    }
    return trans_mat;
  }

  // Store -v_imag if v_imag is not approximately zero. Use this to check if
  // the complex conjugate of an irrep has already been included in the
  // result.
  std::vector<Eigen::VectorXd> real_axes;
  std::vector<Eigen::VectorXd> conj_imag_axes;
  std::vector<bool> found_conj_imag;

  auto add_vector = [&](Eigen::VectorXd const &v) {
    if (row >= trans_mat.rows()) {
      throw std::runtime_error(
          "Error in full_trans_mat: row out of range error");
    }
    trans_mat.block(row, 0, 1, col) = v.transpose();
    row += 1;
  };

  double tol = TOL;
  for (auto const &irrep : irreps) {
    for (Index i = 0; i < irrep.irrep_dim; ++i) {
      Eigen::VectorXd v_imag = irrep.trans_mat.row(i).imag().cast<double>();
      v_imag.normalize();
      Eigen::VectorXd v_real = irrep.trans_mat.row(i).real().cast<double>();
      v_real.normalize();

      double imag_norm = v_imag.norm();

      if (imag_norm > tol) {
        bool found = false;
        for (Index j = 0; j < conj_imag_axes.size(); ++j) {
          if (almost_equal(v_imag, conj_imag_axes[j], tol) &&
              almost_equal(v_real, real_axes[j], tol)) {
            found = true;
            found_conj_imag[j] = true;
            break;
          }
        }
        if (!found) {
          real_axes.push_back(v_real);
          conj_imag_axes.push_back(-v_imag);
          found_conj_imag.push_back(false);
          add_vector(v_real);
          add_vector(v_imag);
        }
      } else {
        add_vector(v_real);
      }
    }
  }

  /// Check that all complex conjugate pairs of irrep vectors have been included
  /// in the result
  for (Index i = 0; i < conj_imag_axes.size(); ++i) {
    if (!found_conj_imag[i]) {
      throw std::runtime_error(
          "Error in full_trans_mat: did not find conjugate pair for all "
          "complex irrep vectors");
    }
  }

  return trans_mat;
}

struct SubspaceIrrepDecomposition {
  SubspaceIrrepDecomposition(Eigen::MatrixXd const &_subspace)
      : subspace(_subspace), complete_decomposition(false) {}

  /// \brief The initial subspace. Must be invariant under the group.
  Eigen::MatrixXd subspace;

  // Work variables:

  /// \brief The kernel of the initial subspace.
  Eigen::MatrixXd initial_kernel;

  // Results:

  /// \brief Subspace that could not be decomposed into irreps.
  Eigen::MatrixXd incomplete_subspace;

  /// \brief Irreps found in the subspace.
  std::vector<IrrepInfo> irreps;

  /// \brief True if all irreps in the subspace were found
  bool complete_decomposition;

  void solve(MatrixRep const &fullspace_rep, GroupIndices const &head_group,
             std::optional<GroupIndicesOrbitSet> const &subgroup_orbits,
             bool allow_complex, std::optional<Log> log);
};

/// \brief Solve for irrep subspaces of a subspace
///
/// \param fullspace_rep Full space matrix representation (rep[0].rows() ==
///     subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param subgroup_orbits The orbits of subgroups to use for symmetrization.
///     If not provided, no symmetrization is performed.
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
/// \param log Optional log to write to
///
/// \returns (irreps, complete_decomposition) where irreps is the vector of
///     IrrepInfo for the irreps found in the subspace, and
///     complete_decomposition is true if all irreps in the subspace were
///     found (i.e. the subspace was fully decomposed into irreps), and false
///     if not (i.e. the subspace was only partially decomposed into irreps,
///     and there is a remaining subspace that is invariant under the group
///     but could not be decomposed into irreps).
void SubspaceIrrepDecomposition::solve(
    MatrixRep const &fullspace_rep, GroupIndices const &head_group,
    std::optional<GroupIndicesOrbitSet> const &subgroup_orbits,
    bool allow_complex, std::optional<Log> log) {
  using namespace IrrepDecompositionImpl;

  // expand the initial subspace into an invariant subspace
  subspace = make_invariant_space(fullspace_rep, head_group, subspace);
  incomplete_subspace = subspace;
  initial_kernel = make_kernel(subspace);

  // work variables:
  Index dim = fullspace_rep[0].rows();

  // 2) Perform irrep_decomposition
  // In some cases the `irrep_decomposition` method does not find all irreps.
  // As long as it finds at least one, this loop will try again in the
  // remaining subspace.
  Index i = 1;
  Index rotation_count = 0;
  while (true) {
    if (log.has_value()) {
      std::stringstream ss;
      ss << "Iteration " << i;
      log->begin<Log::standard>(ss.str());
      log->increase_indent();
      log->indent() << std::endl;

      log->indent() << "Subspace:" << std::endl;
      for (Index j = 0; j < subspace.cols(); ++j) {
        log->indent() << "- " << j << ": "
                      << pretty(subspace.col(j)).transpose() << std::endl;
      }
      log->indent() << std::endl;
    }

    // Irreps are found in a subspace specified via the subspace matrix rep
    // the input `incomplete_subspace` must be orthonormal
    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "Begin subspace matrix representation construction";
      append_time(*log, 1);
    }
    MatrixRep subspace_rep_i =
        make_subspace_rep(fullspace_rep, incomplete_subspace);

    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "DONE";
      append_time(*log, 2);
    }

    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex, log);

    // If no irreps found in the subspace,
    // then we stop with an incomplete decomposition
    if (subspace_irreps_i.size() == 0) {
      if (log.has_value()) {
        log->indent() << "Break: No irreps found in subspace." << std::endl
                      << std::endl;
        log->decrease_indent();
      }

      //

      if (rotation_count >= 10) {
        if (log.has_value()) {
          log->indent() << "Break: Maximum number of iterations reached."
                        << std::endl
                        << std::endl;
          log->decrease_indent();
        }
        break;
      } else {
        // If no irreps were found, try again with a randomly rotated subspace.
        {
          Index n_cols = incomplete_subspace.cols();
          std::random_device rd;
          std::mt19937 gen(rd());
          std::normal_distribution<double> dist(0.0, 1.0);
          Eigen::MatrixXd random_matrix(n_cols, n_cols);
          for (Index r = 0; r < n_cols; ++r) {
            for (Index c = 0; c < n_cols; ++c) {
              random_matrix(r, c) = dist(gen);
            }
          }
          Eigen::HouseholderQR<Eigen::MatrixXd> qr(random_matrix);
          incomplete_subspace =
              incomplete_subspace * Eigen::MatrixXd(qr.householderQ());
        }

        if (log.has_value()) {
          log->indent() << "Randomly rotated subspace (attempt "
                        << rotation_count + 1 << " / 10)" << std::endl;
        }

        ++rotation_count;
        continue;
      }
    }

    // Optionally, symmetrize all the irreps that were found
    if (subgroup_orbits.has_value()) {
      if (log.has_value()) {
        log->begin<Log::standard>("Symmetrization");
        log->increase_indent();
        log->indent() << std::endl;
        log->indent() << "Symmetrization = "
                      << (subgroup_orbits.has_value() ? "yes" : "no")
                      << std::endl
                      << std::endl;
      }

      subspace_irreps_i = symmetrize_irreps(
          subspace_rep_i, head_group, subspace_irreps_i, *subgroup_orbits, log);

      if (log.has_value()) {
        log->indent() << "Symmetrization: DONE" << std::endl << std::endl;
        log->decrease_indent();
      }
    }

    std::vector<IrrepInfo> fullspace_irreps_i =
        make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

    // Save the new fullspace irreps
    for (auto const &irrep : fullspace_irreps_i) {
      irreps.push_back(irrep);
    }

    // Combine the initial kernel and irrep spaces to generate the
    // currently finished subspace
    Eigen::MatrixXd finished_subspace = initial_kernel;

    if (irreps.size()) {
      Eigen::MatrixXd new_subspace =
          full_trans_mat(fullspace_irreps_i, allow_complex).adjoint();
      finished_subspace = extend(
          finished_subspace, full_trans_mat(irreps, allow_complex).adjoint());
    }

    if (finished_subspace.cols() == dim) {
      if (log.has_value()) {
        log->indent() << "Break: Complete decomposition found." << std::endl
                      << std::endl;
        log->decrease_indent();
      }
      complete_decomposition = true;
      break;
    }

    // If not all irreps have been found, try again in remaining space
    incomplete_subspace = make_kernel(finished_subspace);

    if (log.has_value()) {
      log->indent() << "Iteration " << i << ": DONE." << std::endl << std::endl;
      log->decrease_indent();
    }

    // Update iteration count
    ++i;
  }
}

/// IrrepDecomposition constructor
///
/// \param rep Full space matrix representation (rep[0].rows() ==
///     init_subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param init_subspace Input subspace in which irreps are to be found. Will
/// be
///     expanded (column space increased) by application of rep and
///     orthogonalization to form an invariant subspace (i.e. column space
///     dimension is not increased by application of elements in head_group)
/// \param subgroup_orbits The orbits of subgroups to use for symmetrization.
///     If not provided, no symmetrization is performed.
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    std::optional<GroupIndicesOrbitSet> const &subgroup_orbits,
    bool allow_complex, std::optional<Log> _log)
    : init_subspace(init_subspace),
      fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      complete_decomposition(false),
      log(_log) {
  using namespace IrrepDecompositionImpl;

  if (fullspace_rep.size() == 0) {
    std::stringstream msg;
    msg << "Error in IrrepDecomposition: fullspace_rep.size() == 0.";
    throw std::runtime_error(msg.str());
  }

  Index dim = fullspace_rep[0].rows();

  if (log.has_value()) {
    log->begin<Log::standard>("IrrepDecomposition");
    log->indent() << std::endl;
    log->indent() << "Number of elements = " << fullspace_rep.size()
                  << std::endl;
    log->indent() << "Vector space dimension = " << dim << std::endl
                  << std::endl;

    prettyp<Log::debug>(*log, "1. Initial vector space", init_subspace);
    log->indent() << std::endl;

    log->indent() << "Make invariant vector space..." << std::endl;
  }

  // 1) Expand subspace by application of group, and orthonormalization
  subspace = make_invariant_space(fullspace_rep, head_group, init_subspace);
  initial_kernel = make_kernel(subspace);

  if (log.has_value()) {
    log->indent() << "Make invariant vector space: DONE" << std::endl
                  << std::endl;
    log->indent() << "Invariant vector space dimension = " << subspace.cols()
                  << std::endl
                  << std::endl;

    prettyp<Log::debug>(*log, "2. Initial invariant vector space", subspace);
  }

  if (subspace.cols() == 0) {
    std::stringstream msg;
    msg << "Error in IrrepDecomposition: invariant subspace has zero "
           "dimension.";
    throw std::runtime_error(msg.str());
  }

  // 2) Perform irrep_decomposition
  // In some cases the `irrep_decomposition` method does not find all irreps.
  // As long as it finds at least one, this loop will try again in the
  // remaining subspace.
  incomplete_subspace = subspace;
  Index i = 1;
  while (true) {
    if (log.has_value()) {
      log->indent() << "-- Begin iteration " << i;
      append_time(*log, 2);
    }

    // Irreps are found in a subspace specified via the subspace matrix rep
    // the input `incomplete_subspace` must be orthonormal
    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "Begin subspace matrix representation construction";
      append_time(*log, 1);
    }
    MatrixRep subspace_rep_i =
        make_subspace_rep(fullspace_rep, incomplete_subspace);
    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "DONE";
      append_time(*log, 2);
    }

    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex, log);
    if (log.has_value()) {
      print_irreps<Log::debug>(*log, "Irreps, as found", subspace_irreps_i);
    }

    // If no irreps found in the subspace,
    // then we stop with an incomplete decomposition
    if (subspace_irreps_i.size() == 0) {
      // std::stringstream msg;
      // msg << "Error in IrrepDecomposition: failed to find all irreps";
      // throw std::runtime_error(msg.str());
      if (log.has_value()) {
        log->indent() << "Break: No irreps found in subspace." << std::endl
                      << std::endl;
      }

      break;
    }

    // Symmetrize all the irreps that were found
    std::vector<IrrepInfo> fullspace_irreps_i;
    if (log.has_value()) {
      log->begin<Log::standard>("Symmetrization");
      log->indent() << std::endl;
    }

    if (!subgroup_orbits.has_value()) {
      if (log.has_value()) {
        log->indent() << "Symmetrization = no" << std::endl << std::endl;
      }

      fullspace_irreps_i =
          make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

      if (log.has_value()) {
        print_irreps<Log::debug>(*log, "Irreps, full dim", fullspace_irreps_i);
      }
    } else {
      if (log.has_value()) {
        log->indent() << "Symmetrization = yes" << std::endl << std::endl;
      }

      subspace_irreps_i = symmetrize_irreps(
          subspace_rep_i, head_group, subspace_irreps_i, *subgroup_orbits, log);

      if (log.has_value()) {
        log->indent() << "Symmetrization: DONE" << std::endl << std::endl;

        print_irreps<Log::debug>(*log, "Irreps, symmetrized",
                                 subspace_irreps_i);
      }

      fullspace_irreps_i =
          make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

      if (log.has_value()) {
        print_irreps<Log::debug>(*log, "Irreps, symmetrized and full dim",
                                 fullspace_irreps_i);
      }
    }

    // Save the new fullspace irreps
    for (auto const &irrep : fullspace_irreps_i) {
      irreps.push_back(irrep);
    }
    if (log.has_value()) {
      print_irreps<Log::debug>(*log, "Irreps, all found so far",
                               subspace_irreps_i);
    }

    // Combine the initial kernel and irrep spaces to generate the
    // currently finished subspace
    Eigen::MatrixXd finished_subspace = initial_kernel;

    if (irreps.size()) {
      finished_subspace = extend(
          finished_subspace, full_trans_mat(irreps, allow_complex).adjoint());
    }

    if (log.has_value()) {
      Eigen::MatrixXd finished_subspace_i =
          full_trans_mat(fullspace_irreps_i, allow_complex).adjoint();

      prettyp<Log::debug>(*log, "Combined vector space, this step",
                          finished_subspace_i);
    }

    if (log.has_value()) {
      prettyp<Log::debug>(*log, "Combined vector space, so far",
                          finished_subspace);
    }

    if (finished_subspace.cols() == dim) {
      if (log.has_value()) {
        log->indent() << "Break: Complete decomposition found." << std::endl
                      << std::endl;
      }
      complete_decomposition = true;
      break;
    }

    // If not all irreps have been found, try again in remaining space
    incomplete_subspace = make_kernel(finished_subspace);
    if (log.has_value()) {
      log->indent() << " Incomplete decomposition, "
                    << incomplete_subspace.cols() << " / " << subspace.cols()
                    << " dimensions remaining." << std::endl
                    << std::endl;

      prettyp<Log::debug>(*log, "Remaining vector space", incomplete_subspace);
    }

    // Check iteration count
    ++i;
  }

  // 3) Combine to form symmetry adapted subspace
  symmetry_adapted_subspace = full_trans_mat(irreps, allow_complex).adjoint();
  if (log.has_value()) {
    print_irreps<Log::debug>(*log, "3. Irreps, symmetry adapted", irreps);
    prettyp<Log::debug>(*log, "4. Symmetry adapted vector space",
                        symmetry_adapted_subspace);
  }
}

std::set<std::set<Index>> make_disjoint_variable_sets(
    MatrixRep const &fullspace_rep, GroupIndices const &head_group,
    double zero_tol) {
  std::set<std::set<Index>> disjoint_variable_sets;
  std::vector<bool> variable_included(fullspace_rep[0].rows(), false);
  for (Index i_x = 0; i_x < fullspace_rep[0].rows(); ++i_x) {
    if (variable_included[i_x]) {
      continue;
    }

    std::set<Index> variable_set;
    variable_set.insert(i_x);
    variable_included[i_x] = true;
    for (Index j_x = 0; j_x < fullspace_rep[0].rows(); ++j_x) {
      if (variable_included[j_x]) {
        continue;
      }
      for (Index i_g : head_group) {
        if (!almost_zero(fullspace_rep[i_g](j_x, i_x), zero_tol)) {
          variable_set.insert(j_x);
          variable_included[j_x] = true;
          break;
        }
      }
    }
    disjoint_variable_sets.insert(variable_set);
  }
  return disjoint_variable_sets;
}

/// \brief Project a general subspace onto the subspace of variables
///     specified by variable_set and orthonormalize.
///
/// \param variable_set Set of variable indices to project onto
/// \param subspace Subspace to project onto variable set (subspace.rows() ==
///     full space dimension, subspace.cols() == subspace dimension)
/// \param zero_tol Tolerance for determining if the projection is zero.
///
/// \returns Orthonormal basis for the projection of subspace onto the variables
/// specified by variable_set. If the projection is zero, then an empty matrix
/// with shape (subspace.rows(), 0) is returned.
///
Eigen::MatrixXd project_onto_variable_set(Eigen::MatrixXd const &subspace,
                                          std::set<Index> const &variable_set,
                                          double zero_tol) {
  if (variable_set.size() == 0) {
    throw std::runtime_error(
        "Error in project_onto_variable_set: variable_set is empty.");
  }
  if (subspace.cols() == 0) {
    throw std::runtime_error(
        "Error in project_onto_variable_set: subspace has zero columns.");
  }
  if (subspace.rows() == 0) {
    throw std::runtime_error(
        "Error in project_onto_variable_set: subspace has zero rows.");
  }
  Index dim = subspace.rows();

  Eigen::MatrixXd B =
      Eigen::MatrixXd::Zero(variable_set.size(), subspace.cols());
  for (Index i_col = 0; i_col < subspace.cols(); ++i_col) {
    Index idx = 0;
    for (Index j_x : variable_set) {
      B(idx, i_col) = subspace(j_x, i_col);
      ++idx;
    }
  }

  if (B.cwiseAbs().maxCoeff() < zero_tol) {
    return Eigen::MatrixXd::Zero(dim, 0);
  }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(B);
  qr.setThreshold(zero_tol);
  Eigen::MatrixXd Q = qr.householderQ();
  Index rank = qr.rank();
  // This is probably redundant.
  if (rank == 0) {
    return Eigen::MatrixXd::Zero(dim, 0);
  }

  Eigen::MatrixXd Qi = Q.leftCols(rank);
  Qi = standardize_column_vector_signs(Qi, TOL);
  Eigen::MatrixXd subspace_projection = Eigen::MatrixXd::Zero(dim, rank);
  for (Index i_col = 0; i_col < rank; ++i_col) {
    Index idx = 0;
    for (Index j_x : variable_set) {
      subspace_projection(j_x, i_col) = Qi(idx, i_col);
      ++idx;
    }
  }
  return subspace_projection;
};

/// IrrepDecomposition constructor
///
/// \param rep Full space matrix representation (rep[0].rows() ==
///     init_subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param init_subspace Input subspace in which irreps are to be found. Will
/// be
///     expanded (column space increased) by application of rep and
///     orthogonalization to form an invariant subspace (i.e. column space
///     dimension is not increased by application of elements in head_group)
/// \param subgroup_orbits The orbits of subgroups to use for symmetrization.
///     If not provided, no symmetrization is performed.
/// \param class_indices The conjugacy class indices of the group elements,
///     for printing the character table.
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
/// \param flag Select this constructor
/// \param zero_tol Tolerance for determining if variables mix (i.e. if any
///     element of the group representation is nonzero up to this tolerance,
///     then the corresponding variables are considered to mix). This is used
///     to break the input subspace into smaller subspaces that can be solved
///     separately.
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    std::optional<GroupIndicesOrbitSet> const &subgroup_orbits,
    std::optional<std::vector<Index>> const &class_indices, bool allow_complex,
    std::optional<Log> _log, SolveByDisjointVariableSetsFlag /*flag*/,
    double zero_tol)
    : init_subspace(init_subspace),
      fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      complete_decomposition(false),
      log(_log) {
  using namespace IrrepDecompositionImpl;

  if (fullspace_rep.size() == 0) {
    std::stringstream msg;
    msg << "Error in IrrepDecomposition: fullspace_rep.size() == 0.";
    throw std::runtime_error(msg.str());
  }

  Index dim = fullspace_rep[0].rows();

  // Check for variables that do not mix
  std::set<std::set<Index>> disjoint_variable_sets =
      make_disjoint_variable_sets(fullspace_rep, head_group, zero_tol);

  if (log.has_value()) {
    log->begin<Log::standard>("IrrepDecomposition");
    log->indent() << std::endl;
    log->indent() << "Number of elements = " << fullspace_rep.size()
                  << std::endl;
    log->indent() << "Vector space dimension = " << dim << std::endl
                  << std::endl;
    log->indent() << "Number of disjoint variable sets = "
                  << disjoint_variable_sets.size() << std::endl
                  << std::endl;
    log->indent() << "Variable sets:" << std::endl;
    for (const auto &variable_set : disjoint_variable_sets) {
      log->indent() << "- " << SetPrinter(variable_set) << std::endl;
    }
    log->indent() << std::endl;
    log->indent() << "Make invariant vector space..." << std::endl;
  }

  // 1) Expand subspace by application of group, and orthonormalization
  subspace = make_invariant_space(fullspace_rep, head_group, init_subspace);
  initial_kernel = make_kernel(subspace);
  incomplete_subspace = subspace;

  if (log.has_value()) {
    log->indent() << "Make invariant vector space: DONE" << std::endl
                  << std::endl;
    log->indent() << "Invariant vector space dimension = " << subspace.cols()
                  << std::endl
                  << std::endl;
  }

  if (subspace.cols() == 0) {
    std::stringstream msg;
    msg << "Error in IrrepDecomposition: invariant subspace has zero "
           "dimension.";
    throw std::runtime_error(msg.str());
  }

  /// Find irreps for the projection of the initial subspace onto each
  /// variable set subspace
  Index i_variable_set = 1;
  bool running_complete_decomposition = true;
  for (const auto &variable_set : disjoint_variable_sets) {
    // Log variable set
    if (log.has_value()) {
      std::stringstream ss;
      ss << "Variable set " << i_variable_set << " / "
         << disjoint_variable_sets.size();
      log->begin<Log::standard>(ss.str());
      log->indent() << std::endl;
      log->increase_indent();
      log->indent() << "Variable set size = " << variable_set.size()
                    << std::endl;
      log->indent() << "Variable set: " << SetPrinter(variable_set) << std::endl
                    << std::endl;
    }

    Eigen::MatrixXd subspace_i =
        project_onto_variable_set(subspace, variable_set, zero_tol);

    if (subspace_i.cols() == 0) {
      if (log.has_value()) {
        log->indent() << "No overlap with input subspace, skipping."
                      << std::endl
                      << std::endl;
        log->decrease_indent();
      }
      ++i_variable_set;
      continue;
    }

    // // debug
    // if (subspace_i.cols() > B.cols()) {
    //   throw std::runtime_error(
    //       "Error in IrrepDecomposition: subspace_i.cols() > B.cols()");
    // }

    SubspaceIrrepDecomposition x(subspace_i);
    x.solve(fullspace_rep, head_group, subgroup_orbits, allow_complex, log);

    Eigen::MatrixXd symmetry_adapted_subspace_i =
        full_trans_mat(x.irreps, allow_complex).adjoint();

    irreps.insert(irreps.end(), x.irreps.begin(), x.irreps.end());

    if (!x.complete_decomposition) {
      running_complete_decomposition = false;
    }
    ++i_variable_set;

    if (log.has_value()) {
      log->decrease_indent();
    }
  }

  // Finalize results:
  complete_decomposition = running_complete_decomposition;

  // 3) Combine to form symmetry adapted subspace
  Eigen::MatrixXd finished_subspace = initial_kernel;
  if (irreps.size()) {
    finished_subspace = extend(finished_subspace,
                               full_trans_mat(irreps, allow_complex).adjoint());
  }
  symmetry_adapted_subspace = full_trans_mat(irreps, allow_complex).adjoint();
  incomplete_subspace = make_kernel(finished_subspace);

  if (log.has_value()) {
    log->indent() << "IrrepDecomposition: DONE" << std::endl << std::endl;

    log->indent() << "Complete decomposition = " << std::boolalpha
                  << complete_decomposition << std::endl
                  << std::endl;

    log->indent() << "Found " << symmetry_adapted_subspace.cols() << " / "
                  << subspace.cols() << " dimensions." << std::endl
                  << std::endl;

    log->indent() << "Number of irreps found = " << irreps.size() << std::endl
                  << std::endl;
  }
}

}  // namespace irreps
}  // namespace CASM
