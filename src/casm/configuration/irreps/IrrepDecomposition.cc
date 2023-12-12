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
/// \param _cyclic_subgroups Cyclic subgroups of head_group. Cyclic subgropus
///     are those formed by repeated application of a single element. Used for
///     symmetrization of the irrep subspaces.
/// \param _all_subgroups All subgroups of head_group. Used for
///     symmetrization of the irrep subspaces if symmetrization using
///     _cyclic_subgroups fails.
/// \param allow_complex If true, all irreps may be complex-valued, if false,
///     complex irreps are combined to form real representations
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    std::function<GroupIndicesOrbitSet()> make_cyclic_subgroups_f,
    std::function<GroupIndicesOrbitSet()> make_all_subgroups_f,
    bool allow_complex, std::optional<Log> _log)
    : fullspace_rep(_fullspace_rep), head_group(_head_group), log(_log) {
  using namespace IrrepDecompositionImpl;

  if (log.has_value()) {
    log->begin<Log::verbose>("IrrepDecomposition");
    log->indent() << std::endl;
    prettyp<Log::verbose>(*log, "1. Initial vector space", init_subspace);
  }

  Index dim = fullspace_rep[0].rows();

  // 1) Expand subspace by application of group, and orthonormalization
  subspace = make_invariant_space(fullspace_rep, head_group, init_subspace);
  if (log.has_value()) {
    prettyp<Log::verbose>(*log, "2. Initial invariant vector space", subspace);
  }

  // 2) Perform irrep_decomposition
  // In some cases the `irrep_decomposition` method does not find all irreps.
  // As long as it finds at least one, this loop will try again in the remaining
  // subspace.
  Eigen::MatrixXd subspace_i = subspace;
  Eigen::MatrixXd finished_subspace = make_kernel(subspace);
  Index i = 1;
  while (finished_subspace.cols() != dim) {
    if (log.has_value() && log->print()) {
      log->indent() << std::endl;
      log->indent() << "-- Begin step " << i << " --" << std::endl;
    }

    // Irreps are found in a subspace specified via the subspace matrix rep
    MatrixRep subspace_rep_i = make_subspace_rep(fullspace_rep, subspace_i);
    std::vector<IrrepInfo> subspace_irreps_i =
        irrep_decomposition(subspace_rep_i, head_group, allow_complex);
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
    std::vector<IrrepInfo> symmetrized_subspace_irreps_i =
        symmetrize_irreps(subspace_rep_i, head_group, subspace_irreps_i,
                          make_cyclic_subgroups_f, make_all_subgroups_f);
    if (log.has_value()) {
      print_irreps<Log::verbose>(*log, "Irreps, symmetrized",
                                 subspace_irreps_i);
    }

    // Transform the irreps trans_mat to act on vectors in the fullspace
    std::vector<IrrepInfo> symmetrized_fullspace_irreps_i =
        make_fullspace_irreps(symmetrized_subspace_irreps_i, subspace_i);
    if (log.has_value()) {
      print_irreps<Log::verbose>(*log, "Irreps, symmetrized and full dim",
                                 subspace_irreps_i);
    }
    // Save the new fullspace irreps
    for (auto const &irrep : symmetrized_fullspace_irreps_i) {
      irreps.push_back(irrep);
    }
    if (log.has_value()) {
      print_irreps<Log::verbose>(*log, "Irreps, all found so far",
                                 subspace_irreps_i);
    }

    // Combine the irrep spaces and add to finished_subspace
    Eigen::MatrixXd finished_subspace_i =
        full_trans_mat(symmetrized_fullspace_irreps_i).adjoint();
    if (log.has_value()) {
      prettyp<Log::verbose>(*log, "Combined vector space, this step",
                            finished_subspace_i);
    }

    finished_subspace = extend(finished_subspace, finished_subspace_i);
    if (log.has_value()) {
      prettyp<Log::verbose>(*log, "Combined vector space, so far",
                            finished_subspace_i);
    }

    // If not all irreps have been found, try again in remaining space
    subspace_i = make_kernel(finished_subspace);
    if (log.has_value()) {
      prettyp<Log::verbose>(*log, "Remaining vector space", subspace_i);
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

}  // namespace irreps
}  // namespace CASM
