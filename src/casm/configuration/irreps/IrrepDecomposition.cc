#include "casm/configuration/irreps/IrrepDecomposition.hh"

#include <iostream>

#include "casm/configuration/irreps/IrrepDecompositionImpl.hh"
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
/// \param make_cyclic_subgroups_f Function that returns cyclic subgroups of
///     head_group. Cyclic subgropus are those formed by repeated application
///     of a single element. Used for symmetrization of the irrep subspaces if
///     symmetrization == "fast".
/// \param make_all_subgroups_f Function that returns all subgroups of
///     head_group. Used for symmetrization of the irrep subspaces if
///     symmetrization == "complete".
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
/// \param symmetrization Type of symmetrization to perform on irrep subspaces
///     Options:
///     - "none": Leave the irreducible subspace bases as initially found,
///       reducing computation time.
///     - "fast": Symmetrize the irreducible subspace bases to align along
///       high-symmetry directions using cyclic subgroups. This may not be a
///       complete symmetrization, but is generally fast.
///     - "complete": Symmetrize the irreducible subspace bases to align
///       along high-symmetry directions using all subgroups. For large
///       spaces, finding all subgroups is slow.
/// \param max_iter Irrep decomposition is not guaranteed to be complete for
///     a particular choice of `init_subspace`. When it is incomplete, the
///     irrep decomposition is repeated on the remaining subspace. This
///     continues for up to a maximum number of iterations, specified by
///     `max_iter`.
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    std::function<GroupIndicesOrbitSet()> make_cyclic_subgroups_f,
    std::function<GroupIndicesOrbitSet()> make_all_subgroups_f,
    bool allow_complex, std::string _symmetrization, Index max_iter,
    std::optional<Log> _log)
    : init_subspace(init_subspace),
      symmetrization(_symmetrization),
      fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      complete_decomposition(false),
      log(_log) {
  using namespace IrrepDecompositionImpl;

  std::vector<std::string> symm_options = {"none", "fast", "complete"};
  if (std::find(symm_options.begin(), symm_options.end(), symmetrization) ==
      symm_options.end()) {
    std::stringstream msg;
    msg << "Error in IrrepDecomposition: invalid symmetrization option '"
        << symmetrization << "'. Valid options are: ";
    for (auto const &opt : symm_options) {
      msg << "'" << opt << "' ";
    }
    throw std::runtime_error(msg.str());
  }

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

    prettyp<Log::verbose>(*log, "1. Initial vector space", init_subspace);
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

    prettyp<Log::verbose>(*log, "2. Initial invariant vector space", subspace);
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
      log->indent() << "-- Begin iteration " << i << " --" << std::endl
                    << std::endl;
    }

    // Irreps are found in a subspace specified via the subspace matrix rep
    MatrixRep subspace_rep_i =
        make_subspace_rep(fullspace_rep, incomplete_subspace);
    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex, log);
    if (log.has_value()) {
      print_irreps<Log::verbose>(*log, "Irreps, as found", subspace_irreps_i);
    }

    // If not irreps found in the subspace, this method has failed
    // If the irreps do not span the whole subspace, we'll try again
    if (subspace_irreps_i.size() == 0) {
      std::stringstream msg;
      msg << "Error in IrrepDecomposition: failed to find all irreps";
      throw std::runtime_error(msg.str());
    }

    // Symmetrize all the irreps that were found
    std::vector<IrrepInfo> fullspace_irreps_i;
    if (log.has_value()) {
      log->increase_indent();
      log->begin<Log::standard>("Symmetrization");
      log->indent() << std::endl;
    }

    if (symmetrization == "none") {
      if (log.has_value()) {
        log->indent() << "Symmetrization = None" << std::endl << std::endl;
      }

      fullspace_irreps_i =
          make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

      if (log.has_value()) {
        print_irreps<Log::verbose>(*log, "Irreps, full dim",
                                   fullspace_irreps_i);
      }
    } else if (symmetrization == "fast") {
      if (log.has_value()) {
        log->indent() << "Symmetrization = fast" << std::endl << std::endl;
        log->indent() << "Begin fast symmetrization..." << std::endl
                      << std::endl;
      }

      subspace_irreps_i =
          symmetrize_irreps(subspace_rep_i, head_group, subspace_irreps_i,
                            make_cyclic_subgroups_f, log);

      if (log.has_value()) {
        log->indent() << "Fast symmetrization: DONE" << std::endl << std::endl;

        print_irreps<Log::verbose>(*log, "Irreps, symmetrized",
                                   subspace_irreps_i);
      }

      fullspace_irreps_i =
          make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

      if (log.has_value()) {
        print_irreps<Log::verbose>(*log, "Irreps, symmetrized and full dim",
                                   fullspace_irreps_i);
      }
    } else if (symmetrization == "complete") {
      if (log.has_value()) {
        log->indent() << "Symmetrization = complete" << std::endl << std::endl;
        log->indent() << "Begin complete symmetrization..." << std::endl
                      << std::endl;
      }

      subspace_irreps_i =
          symmetrize_irreps(subspace_rep_i, head_group, subspace_irreps_i,
                            make_all_subgroups_f, log);

      if (log.has_value()) {
        log->indent() << "Complete symmetrization: DONE" << std::endl
                      << std::endl;

        print_irreps<Log::verbose>(*log, "Irreps, symmetrized",
                                   subspace_irreps_i);
      }

      fullspace_irreps_i =
          make_fullspace_irreps(subspace_irreps_i, incomplete_subspace);

      if (log.has_value()) {
        print_irreps<Log::verbose>(*log, "Irreps, symmetrized and full dim",
                                   fullspace_irreps_i);
      }
    }

    if (log.has_value()) {
      log->end_section();  // symmetrization
      log->decrease_indent();
    }

    // Save the new fullspace irreps
    for (auto const &irrep : fullspace_irreps_i) {
      irreps.push_back(irrep);
    }
    if (log.has_value()) {
      print_irreps<Log::verbose>(*log, "Irreps, all found so far",
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

      prettyp<Log::verbose>(*log, "Combined vector space, this step",
                            finished_subspace_i);
    }

    if (log.has_value()) {
      prettyp<Log::verbose>(*log, "Combined vector space, so far",
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

      prettyp<Log::verbose>(*log, "Remaining vector space",
                            incomplete_subspace);
    }

    // Check iteration count
    ++i;
    if (i > max_iter) {
      if (log.has_value()) {
        log->indent() << "Break: Maximum number of iterations (" << max_iter
                      << ") reached." << std::endl
                      << std::endl;
      }
      break;
    }
  }

  // 3) Combine to form symmetry adapted subspace
  symmetry_adapted_subspace = full_trans_mat(irreps).adjoint();
  if (log.has_value()) {
    print_irreps<Log::verbose>(*log, "3. Irreps, symmetry adapted", irreps);
    prettyp<Log::verbose>(*log, "4. Symmetry adapted vector space",
                          symmetry_adapted_subspace);
  }
}

/// IrrepDecomposition constructor (from existing decomposition)
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &_init_subspace, Eigen::MatrixXd const &_subspace,
    std::vector<IrrepInfo> const &_irreps, bool _complete_decomposition,
    Eigen::MatrixXd const &_incomplete_subspace, std::optional<Log> _log)
    : fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      init_subspace(_init_subspace),
      subspace(_subspace),
      irreps(_irreps),
      complete_decomposition(_complete_decomposition),
      incomplete_subspace(_incomplete_subspace),
      log(_log) {
  using namespace IrrepDecompositionImpl;

  // construct initial kernel:
  initial_kernel = make_kernel(subspace);

  // construct symmetry adapted subspace:
  symmetry_adapted_subspace = full_trans_mat(irreps).adjoint();
}

void IrrepDecomposition::symmetrize_all_irreps(std::string symmetrization) {}

void IrrepDecomposition::symmetrize_irrep(Index i, std::string symmetrization) {
}

}  // namespace irreps
}  // namespace CASM
