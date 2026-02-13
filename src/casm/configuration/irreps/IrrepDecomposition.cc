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

/// \brief Combine individual possibly complex irrep trans_mat to form one
/// larger real trans_mat
///
/// For allow_complex=false, simply concatenate the real parts of the individual
/// trans_mats. For allow_complex=true, add the original Re and Im parts of the
/// individual trans_mats, but avoid duplication from complex conjugate pairs
/// by checking if the real subspace of a complex irrep is already covered by
/// previously added irreps.
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

  // For allow_complex=true:
  // Complex conjugate irrep pairs span the same real subspace. For each
  // complex irrep, check if its real subspace is already covered; if not,
  // add the original Re and Im parts of its rows (preserving symmetrized
  // axes). Skip conjugate partners whose subspace is already covered.
  using IrrepDecompositionImpl::extend;
  using IrrepDecompositionImpl::is_extended_by;

  // Track covered complex subspace (columns are orthonormal basis vectors)
  Eigen::MatrixXcd covered(col, 0);

  for (auto const &irrep : irreps) {
    if (!irrep.complex) {
      // Real irrep: add trans_mat rows directly
      trans_mat.block(row, 0, irrep.irrep_dim, irrep.vector_dim) =
          irrep.trans_mat.real();
      row += irrep.irrep_dim;
    } else {
      // Complex irrep: collect Re and Im parts of all rows as column
      // vectors, then orthogonalize to build a basis for the subspace check
      Eigen::MatrixXd parts(col, 2 * irrep.irrep_dim);
      for (Index i = 0; i < irrep.irrep_dim; ++i) {
        parts.col(2 * i) =
            Eigen::VectorXd(irrep.trans_mat.row(i).real().transpose());
        parts.col(2 * i + 1) =
            Eigen::VectorXd(irrep.trans_mat.row(i).imag().transpose());
      }

      // Orthogonalize for subspace check only
      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> colqr(parts);
      colqr.setThreshold(TOL);
      Index rank = colqr.rank();

      Eigen::HouseholderQR<Eigen::MatrixXd> qr(parts);
      Eigen::MatrixXd Q = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);

      // Check if this real subspace extends the covered subspace.
      // Conjugate partners span the same real subspace, so the second
      // one encountered will not extend and will be skipped.
      Eigen::MatrixXcd Q_complex = Q.cast<std::complex<double>>();
      if (is_extended_by(covered, Q_complex)) {
        // Add Re and Im parts of trans_mat rows, normalized to unit length.
        for (Index i = 0; i < irrep.irrep_dim; ++i) {
          if (row + 1 >= trans_mat.rows()) {
            throw std::runtime_error(
                "Error in full_trans_mat: row out of range error");
          }
          Eigen::RowVectorXd re_row =
              irrep.trans_mat.row(i).real().template cast<double>();
          trans_mat.block(row, 0, 1, col) = re_row.normalized();
          row += 1;
          Eigen::RowVectorXd im_row =
              irrep.trans_mat.row(i).imag().template cast<double>();
          trans_mat.block(row, 0, 1, col) = im_row.normalized();
          row += 1;
        }
        covered = extend(covered, Q_complex);
      }
      // else: conjugate partner already covered, skip
    }
  }

  if (row != trans_mat.rows()) {
    std::stringstream msg;
    msg << "Error in full_trans_mat: expected " << trans_mat.rows()
        << " rows but got " << row;
    throw std::runtime_error(msg.str());
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
             bool allow_complex, std::optional<Log> log,
             CommuterMethod method = CommuterMethod::deterministic);
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
    bool allow_complex, std::optional<Log> log, CommuterMethod method) {
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
  Index iteration_index = 1;
  Index rotation_count = 0;
  while (true) {
    if (log.has_value()) {
      std::stringstream ss;
      ss << "Iteration " << iteration_index;
      log->begin<Log::standard>(ss.str());
      log->increase_indent();
      log->indent() << std::endl;

      // log->indent() << "Incomplete subspace:" << std::endl;
      // for (Index j = 0; j < incomplete_subspace.cols(); ++j) {
      //   log->indent() << "- " << j << ": "
      //                 << pretty(incomplete_subspace.col(j)).transpose()
      //                 << std::endl;
      // }
      // log->indent() << std::endl;
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

    // iteration_index start at 1 ->
    // start with real seed for iteration 1, complex seed for iteration 2, etc.
    bool start_with_real_seed = (iteration_index % 2 == 1);
    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex, log,
                            method, start_with_real_seed);

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
      log->indent() << "Iteration " << iteration_index << ": DONE." << std::endl
                    << std::endl;
      log->decrease_indent();
    }

    // Update iteration count
    ++iteration_index;
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
    bool allow_complex, std::optional<Log> _log, CommuterMethod method)
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
  Index iteration_index = 1;
  while (true) {
    if (log.has_value()) {
      log->indent() << "-- Begin iteration " << iteration_index;
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

    bool start_with_real_seed = (iteration_index % 2 == 1);
    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex, log,
                            method, start_with_real_seed);
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
    ++iteration_index;
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
    double zero_tol, CommuterMethod method)
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

  // Check for columns that do not mix
  MatrixRep subspace_rep = make_subspace_rep(fullspace_rep, subspace);
  std::set<std::set<Index>> disjoint_column_sets =
      make_disjoint_variable_sets(subspace_rep, head_group, zero_tol);

  /// Find irreps for the projection of the initial subspace onto each
  /// variable set subspace
  Index i_variable_set = 1;
  bool running_complete_decomposition = true;
  for (const auto &column_set : disjoint_column_sets) {
    // Log column set
    if (log.has_value()) {
      std::stringstream ss;
      ss << "Column set " << i_variable_set << " / "
         << disjoint_column_sets.size();
      log->begin<Log::standard>(ss.str());
      log->indent() << std::endl;
      log->increase_indent();
      log->indent() << "Column set size = " << column_set.size() << std::endl;
      log->indent() << "Column set: " << SetPrinter(column_set) << std::endl
                    << std::endl;
    }

    // Create subspace_i from columns in column_set
    Eigen::MatrixXd subspace_i =
        Eigen::MatrixXd::Zero(subspace.rows(), column_set.size());
    Index i_col = 0;
    for (Index j : column_set) {
      subspace_i.col(i_col) = subspace.col(j);
      ++i_col;
    }

    if (subspace_i.cols() == 0) {
      if (log.has_value()) {
        log->indent() << "No overlap with input subspace, skipping."
                      << std::endl
                      << std::endl;
        log->decrease_indent();
      }
      ++i_variable_set;
      continue;
    } else {
      if (log.has_value()) {
        log->indent() << "Overlap with input subspace, continue..." << std::endl
                      << std::endl;
      }
    }

    SubspaceIrrepDecomposition x(subspace_i);
    x.solve(fullspace_rep, head_group, subgroup_orbits, allow_complex, log,
            method);

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
    Eigen::MatrixXd irreps_subspace =
        full_trans_mat(irreps, allow_complex).adjoint();
    finished_subspace = extend(finished_subspace, irreps_subspace);
  }
  symmetry_adapted_subspace = full_trans_mat(irreps, allow_complex).adjoint();

  if (finished_subspace.cols() > finished_subspace.rows()) {
    throw std::runtime_error(
        "Error in IrrepDecomposition: finished_subspace has more columns than "
        "rows for unknown reason.");
  }
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
