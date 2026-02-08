#include "casm/configuration/irreps/IrrepDecomposition.hh"

#include "casm/configuration/irreps/IrrepDecompositionImpl.hh"
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
Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps) {
  Index row = 0;
  Index col = 0;
  for (auto const &irrep : irreps) {
    col = irrep.vector_dim;
    row += irrep.irrep_dim;
  }
  Eigen::MatrixXd trans_mat(row, col);
  row = 0;
  for (auto const &irrep : irreps) {
    trans_mat.block(row, 0, irrep.irrep_dim, irrep.vector_dim) =
        irrep.trans_mat.real();
    row += irrep.irrep_dim;
  }
  return trans_mat;
}

/// IrrepDecomposition constructor
///
/// \param rep Full space matrix representation (rep[0].rows() ==
///     init_subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param init_subspace Input subspace in which irreps are to be found. Will be
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
      finished_subspace =
          extend(finished_subspace, full_trans_mat(irreps).adjoint());
    }

    if (log.has_value()) {
      Eigen::MatrixXd finished_subspace_i =
          full_trans_mat(fullspace_irreps_i).adjoint();

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
                    << incomplete_subspace.cols() << " / " << dim
                    << " dimensions remaining." << std::endl
                    << std::endl;

      prettyp<Log::debug>(*log, "Remaining vector space", incomplete_subspace);
    }

    // Check iteration count
    ++i;
  }

  // 3) Combine to form symmetry adapted subspace
  symmetry_adapted_subspace = full_trans_mat(irreps).adjoint();
  if (log.has_value()) {
    print_irreps<Log::debug>(*log, "3. Irreps, symmetry adapted", irreps);
    prettyp<Log::debug>(*log, "4. Symmetry adapted vector space",
                        symmetry_adapted_subspace);
  }
}

}  // namespace irreps
}  // namespace CASM
