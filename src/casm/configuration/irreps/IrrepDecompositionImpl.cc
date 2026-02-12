#include "casm/configuration/irreps/IrrepDecompositionImpl.hh"

#include <iomanip>
#include <iostream>
#include <random>

#include "casm/configuration/irreps/Symmetrizer.hh"
#include "casm/configuration/irreps/misc.hh"
#include "casm/configuration/irreps/to_real.hh"
#include "casm/global/threads.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

namespace irreps {

namespace IrrepDecompositionImpl {

// note: there are a number of floating point comparisons, currently all set
// to use CASM::TOL as the tolerance. If necessary, they could be tuned here or
// via function parameters. To find them search for TOL or "almost".

CommuterParamsCounter::CommuterParamsCounter() : m_valid(false) {}

void CommuterParamsCounter::reset(Eigen::MatrixXcd const &kernel) {
  m_valid = true;
  m_max_cols = kernel.cols();
  m_phase_index = 0;

  kernel_column_i = 0;
  kernel_column_j = m_max_cols - 1;
  phase = std::complex<double>(1., 0.);
}

bool CommuterParamsCounter::valid() const { return m_valid; }

bool CommuterParamsCounter::increment() {
  if (!m_valid) {
    return false;
  }

  // inner loop is column j
  --kernel_column_j;

  // next inner loop is column i
  if (kernel_column_j < kernel_column_i) {
    ++kernel_column_i;
    kernel_column_j = m_max_cols - 1;
  }

  // outer loop is phase: 1, i
  if (kernel_column_i == m_max_cols) {
    ++m_phase_index;
    kernel_column_i = 0;
    kernel_column_j = m_max_cols - 1;

    // once m_phase_index gets to 2, we've finished all possibilities
    if (m_phase_index == 2) {
      m_valid = false;
      return false;
    } else if (m_phase_index == 1) {
      phase = std::complex<double>(0., 1.);
    }
  }

  // skip i==j when phase==i
  if (kernel_column_i == kernel_column_j && m_phase_index == 1) {
    return increment();
  }

  return true;
}

Eigen::MatrixXcd make_commuter(CommuterParamsCounter const &params,
                               MatrixRep const &rep,
                               GroupIndices const &head_group,
                               Eigen::MatrixXcd const &kernel) {
  Index dim = rep[0].rows();

  // construct outer product space of kernel columns i and j
  // commuters are constructed to be self-adjoint,
  // which assures eigenvalues are real
  auto const &col_i = kernel.col(params.kernel_column_i);
  auto const &col_j = kernel.col(params.kernel_column_j);
  auto const &phase = params.phase;
  Eigen::MatrixXcd M_init = phase * col_i * col_j.adjoint() +
                            std::conj(phase) * col_j * col_i.adjoint();

  // If there are no group elements, return zero matrix
  if (head_group.size() == 0) {
    return complex_Zero(dim, dim);
  }

  Eigen::MatrixXcd M = complex_Zero(dim, dim);

  // Convert head_group into an indexable vector for partitioning
  std::vector<Index> head_vec;
  head_vec.reserve(head_group.size());
  for (Index idx : head_group) head_vec.push_back(idx);

  // Determine number of threads to use
  unsigned int hw_conc = std::thread::hardware_concurrency();
  Index max_threads = hw_conc == 0 ? 1 : static_cast<Index>(hw_conc);
  Index n_threads = std::min<Index>(static_cast<Index>(head_vec.size()),
                                    std::max<Index>(1, max_threads));

  // If only one thread, do the simple serial loop for minimal overhead
  if (n_threads == 1) {
    for (Index element_index : head_group) {
      M.noalias() +=
          rep[element_index] * M_init * rep[element_index].transpose();
    }
    return M;
  }

  // Prepare per-thread local accumulators
  std::vector<Eigen::MatrixXcd> local_Ms(n_threads);
  for (Index t = 0; t < n_threads; ++t) local_Ms[t] = complex_Zero(dim, dim);

  // Partition head_vec into contiguous chunks
  Index total = static_cast<Index>(head_vec.size());
  Index chunk = (total + n_threads - 1) / n_threads;

  // Launch threads
  std::vector<std::thread> threads;
  threads.reserve(n_threads);
  for (Index t = 0; t < n_threads; ++t) {
    Index start = t * chunk;
    Index end = std::min(start + chunk, total);

    threads.emplace_back([&, start, end, t]() {
      Eigen::MatrixXcd &local = local_Ms[t];
      for (Index idx = start; idx < end; ++idx) {
        Index element_index = head_vec[idx];
        local.noalias() +=
            rep[element_index] * M_init * rep[element_index].transpose();
      }
    });
  }

  // Join threads
  for (auto &th : threads) {
    if (th.joinable()) th.join();
  }

  // Sum local accumulators into final matrix
  for (Index t = 0; t < n_threads; ++t) {
    M.noalias() += local_Ms[t];
  }

  return M;
}

Eigen::MatrixXcd make_random_commuter(MatrixRep const &rep,
                                      GroupIndices const &head_group,
                                      Eigen::MatrixXcd const &kernel,
                                      std::mt19937 &gen,
                                      bool use_complex_seed) {
  Index dim = rep[0].rows();
  Index k = kernel.cols();

  std::normal_distribution<double> dist(0.0, 1.0);
  Eigen::MatrixXcd M_init;

  if (use_complex_seed) {
    // Generate random complex Hermitian matrix of size (k x k).
    // A complex Hermitian seed can distinguish complex conjugate irrep pairs
    // that share eigenvalues under a real symmetric commuter.
    Eigen::MatrixXcd A(k, k);
    for (Index r = 0; r < k; ++r) {
      for (Index c = 0; c < k; ++c) {
        A(r, c) = std::complex<double>(0.0, dist(gen));
      }
    }
    Eigen::MatrixXcd H = A + A.adjoint();
    M_init = kernel * H * kernel.adjoint();
  } else {
    // Generate random real symmetric matrix of size (k x k).
    // Using a real symmetric seed (rather than complex Hermitian) ensures that
    // the resulting commuter is real symmetric when the kernel and rep are
    // real, which keeps eigenvectors real for real irreps.
    Eigen::MatrixXd A_real(k, k);
    for (Index r = 0; r < k; ++r) {
      for (Index c = 0; c < k; ++c) {
        A_real(r, c) = dist(gen);
      }
    }
    Eigen::MatrixXcd H =
        (A_real + A_real.transpose()).cast<std::complex<double>>();
    M_init = kernel * H * kernel.adjoint();
  }

  // If there are no group elements, return zero matrix
  if (head_group.size() == 0) {
    return complex_Zero(dim, dim);
  }

  Eigen::MatrixXcd M = complex_Zero(dim, dim);

  // Convert head_group into an indexable vector for partitioning
  std::vector<Index> head_vec;
  head_vec.reserve(head_group.size());
  for (Index idx : head_group) head_vec.push_back(idx);

  // Determine number of threads to use
  unsigned int hw_conc = std::thread::hardware_concurrency();
  Index max_threads = hw_conc == 0 ? 1 : static_cast<Index>(hw_conc);
  Index n_threads = std::min<Index>(static_cast<Index>(head_vec.size()),
                                    std::max<Index>(1, max_threads));

  // If only one thread, do the simple serial loop for minimal overhead
  if (n_threads == 1) {
    for (Index element_index : head_group) {
      M.noalias() +=
          rep[element_index] * M_init * rep[element_index].transpose();
    }
    return M;
  }

  // Prepare per-thread local accumulators
  std::vector<Eigen::MatrixXcd> local_Ms(n_threads);
  for (Index t = 0; t < n_threads; ++t) local_Ms[t] = complex_Zero(dim, dim);

  // Partition head_vec into contiguous chunks
  Index total = static_cast<Index>(head_vec.size());
  Index chunk = (total + n_threads - 1) / n_threads;

  // Launch threads
  std::vector<std::thread> threads;
  threads.reserve(n_threads);
  for (Index t = 0; t < n_threads; ++t) {
    Index start = t * chunk;
    Index end = std::min(start + chunk, total);

    threads.emplace_back([&, start, end, t]() {
      Eigen::MatrixXcd &local = local_Ms[t];
      for (Index idx = start; idx < end; ++idx) {
        Index element_index = head_vec[idx];
        local.noalias() +=
            rep[element_index] * M_init * rep[element_index].transpose();
      }
    });
  }

  // Join threads
  for (auto &th : threads) {
    if (th.joinable()) th.join();
  }

  // Sum local accumulators into final matrix
  for (Index t = 0; t < n_threads; ++t) {
    M.noalias() += local_Ms[t];
  }

  return M;
}

Eigen::MatrixXcd make_kernel(Eigen::MatrixXcd const &subspace) {
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr;
  qr.compute(subspace);
  return Eigen::MatrixXcd(qr.householderQ())
      .rightCols(subspace.rows() - subspace.cols());
}
Eigen::MatrixXd make_kernel(Eigen::MatrixXd const &subspace) {
  Eigen::HouseholderQR<Eigen::MatrixXd> qr;
  qr.compute(subspace);
  return Eigen::MatrixXd(qr.householderQ())
      .rightCols(subspace.rows() - subspace.cols());
}

Index find_end_of_equal_eigenvalues(Index begin,
                                    Eigen::VectorXd const &eigenvalues) {
  Index end = begin + 1;
  while (end < eigenvalues.size() &&
         almost_equal(eigenvalues(begin), eigenvalues(end), TOL)) {
    end++;
  }
  return end;
}

/// Make an irreducible subspace
///
/// Makes an orthogonalized irreducible subspace, from (K * V), where
/// - K is the kernel matrix,
/// - V is the eigenvector matrix of ( K.adjoint() * M_new * K )
/// - M_new is the new non-zero commuter matrix
///
/// \param KV_matrix Matrix (K * V)
/// \param begin, end Range of equal eigenvalues
/// \param allow_complex If true, allow subspace with complex basis vectors. If
///     false, will make a pseudo irrep subspace that combines two complex
///     irreps. In this case the irrep is reducible, but this is the most-
///     reduced representation that has real basis vectors.
///
Eigen::MatrixXcd make_irrep_subspace(Eigen::MatrixXcd const &KV_matrix,
                                     Index begin, Index end,
                                     bool allow_complex) {
  Index dim = KV_matrix.rows();

  // get columns of (K * V) corresponding to equal eigenvalues
  Eigen::MatrixXcd X = (KV_matrix).block(0, begin, dim, end - begin);
  Eigen::MatrixXcd subspace_init;
  if (allow_complex) {
    subspace_init = X;
  } else {
    subspace_init = Eigen::MatrixXcd::Zero(dim, 2 * X.cols());
    subspace_init.leftCols(X.cols()) =
        sqrt(2.0) * X.real().cast<std::complex<double>>();
    subspace_init.rightCols(X.cols()) =
        sqrt(2.0) * X.imag().cast<std::complex<double>>();
  }

  // QR decomposition

  // "it seems stupid", but Eigen::HouseholderQR is not rank revealing, and
  // Eigen::ColPivHouseholderQR permutes columns of Q
  Eigen::HouseholderQR<Eigen::MatrixXcd> qr;
  Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> colqr;
  colqr.setThreshold(TOL);

  qr.compute(subspace_init);
  colqr.compute(subspace_init);
  Eigen::MatrixXcd Q = qr.householderQ();
  Eigen::MatrixXcd irrep_subspace = Q.leftCols(colqr.rank());
  return irrep_subspace;
}

/// Calculate character for all matrices in rep
Eigen::VectorXcd make_characters(std::vector<Eigen::MatrixXcd> const &rep) {
  Eigen::VectorXcd characters(rep.size());

  Index element_index = 0;
  for (Eigen::MatrixXcd const &matrix : rep) {
    characters(element_index) = matrix.trace();
    ++element_index;
  }
  return characters;
}

/// Calculate character for all matrices in rep
Eigen::VectorXd make_characters(std::vector<Eigen::MatrixXd> const &rep) {
  Eigen::VectorXd characters(rep.size());

  Index element_index = 0;
  for (Eigen::MatrixXd const &matrix : rep) {
    characters(element_index) = matrix.trace();
    ++element_index;
  }
  return characters;
}

/// Check if approximately zero outside block along diagonal
///
/// Only checks columns and rows in range [begin, end)
bool make_is_block_diagonal(std::vector<Eigen::MatrixXcd> const &rep,
                            Index begin, Index end, double tol) {
  Index n = rep.size();
  std::atomic<bool> is_block_diagonal(true);
  Index len = end - begin;

  auto worker = [&](Index istart, Index iend, Index thread_id) {
    for (Index element_index = istart; element_index < iend; ++element_index) {
      Eigen::MatrixXcd const &matrix = rep[element_index];
      // left
      if (begin != 0) {
        if (!matrix.block(begin, 0, len, begin).isZero(tol)) {
          is_block_diagonal.store(false);
          return;
        }
      }
      // right
      if (end != matrix.cols()) {
        if (!matrix.block(begin, end, len, matrix.cols() - end).isZero(tol)) {
          is_block_diagonal.store(false);
          return;
        }
      }
      // top
      if (begin != 0) {
        if (!matrix.block(0, begin, begin, len).isZero(tol)) {
          is_block_diagonal.store(false);
          return;
        }
      }
      // bottom
      if (end != matrix.rows()) {
        if (!matrix.block(end, begin, matrix.rows() - end, len).isZero(tol)) {
          is_block_diagonal.store(false);
          return;
        }
      }
      if (!is_block_diagonal.load()) {
        return;
      }
    }
  };

  threaded_run(n, worker);

  bool result = is_block_diagonal.load();

  if (!result) {
    throw std::runtime_error(
        "!!! TEST = Representation is not block diagonal !!!");
  }

  return result;
}

/// Find characters for block in range [begin, end)
Eigen::VectorXcd make_irrep_characters(std::vector<Eigen::MatrixXcd> const &rep,
                                       Index begin, Index end) {
  Index n = rep.size();
  Eigen::VectorXcd characters(n);
  Index len = end - begin;

  Index element_index = 0;
  for (Eigen::MatrixXcd const &matrix : rep) {
    characters(element_index) = matrix.block(begin, begin, len, len).trace();
    ++element_index;
  }

  return characters;
}

double make_squared_norm(Eigen::VectorXcd const &characters) {
  double squared_norm = 0.0;
  for (Index i = 0; i < characters.size(); ++i) {
    squared_norm += std::norm(characters(i));
  }
  return squared_norm;
}

double make_squared_norm(Eigen::VectorXd const &characters) {
  double squared_norm = 0.0;
  for (Index i = 0; i < characters.size(); ++i) {
    squared_norm += characters(i) * characters(i);
  }
  return squared_norm;
}

std::complex<double> frobenius_product(Eigen::MatrixXcd const &matrix) {
  return (matrix.array().conjugate() * matrix.array()).sum();
}

Eigen::MatrixXcd normalize_commuter(Eigen::MatrixXcd const &commuter) {
  return commuter / sqrt(frobenius_product(commuter).real());
}

/// Return true if space_A is extended by space_B
bool is_extended_by(Eigen::MatrixXcd const &space_A,
                    Eigen::MatrixXcd const &space_B) {
  return almost_zero((space_B.adjoint() * space_A).norm(), TOL);
}

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXcd extend(Eigen::MatrixXcd const &space_A,
                        Eigen::MatrixXcd const &space_B) {
  Eigen::MatrixXcd result(space_A.rows(), space_A.cols() + space_B.cols());
  result.leftCols(space_A.cols()) = space_A;
  result.rightCols(space_B.cols()) = space_B;
  return result;
}

/// Return matrix combining columns of space_A and space_B
Eigen::MatrixXd extend(Eigen::MatrixXd const &space_A,
                       Eigen::MatrixXd const &space_B) {
  Eigen::MatrixXd result(space_A.rows(), space_A.cols() + space_B.cols());
  result.leftCols(space_A.cols()) = space_A;
  result.rightCols(space_B.cols()) = space_B;
  return result;
}

Index get_total_dim(std::set<PossibleIrrep> const &irreps) {
  Index total_dim = 0;
  for (PossibleIrrep const &irrep : irreps) {
    total_dim += irrep.irrep_dim;
  }
  return total_dim;
}

/// \brief Constructor
///
/// \param eigenvalues Eigenvalues of (K.adjoint() * M_new * K)
/// \param KV_matrix K * V, where V is the eigenvector matrix of
///     (K.adjoint() * M_new * K)
/// \param transformed_rep Transformed representation matrices
/// \param _is_block_diagonal True if representation matrices are block diagonal
/// \param _head_group_size Size of head group
/// \param allow_complex If true, allow subspace with complex basis vectors. If
///     false, will make a pseudo irrep subspace that combines two complex
///     irreps. In this case the irrep is reducible, but this is the most-
///     reduced representation that has real basis vectors.
/// \param _begin, _end Range of columns with equal eigenvalues
/// ///
PossibleIrrep::PossibleIrrep(
    Eigen::VectorXd const &eigenvalues, Eigen::MatrixXcd const &KV_matrix,
    std::vector<Eigen::MatrixXcd> const &transformed_rep,
    bool _is_block_diagonal, Index _head_group_size, bool allow_complex,
    Index _begin, Index _end)
    : head_group_size(_head_group_size),
      begin(_begin),
      end(_end),
      irrep_dim(end - begin),
      is_block_diagonal(_is_block_diagonal) {
  // Log &log = CASM::log();
  // log.indent() << "b";
  // append_time(log, 1);
  characters = make_irrep_characters(transformed_rep, begin, end);

  // log.indent() << "c";
  // append_time(log, 1);
  characters_squared_norm = make_squared_norm(characters);

  bool char_squared_norm_is_head_group_size =
      almost_equal(characters_squared_norm, double(head_group_size), TOL);

  if (!char_squared_norm_is_head_group_size) {
    is_irrep = false;
    return;
  }

  // is_block_diagonal = make_is_block_diagonal(transformed_rep, begin, end,
  // TOL);

  is_irrep = is_block_diagonal && char_squared_norm_is_head_group_size;

  if (!is_irrep) {
    return;
  }

  // log.indent() << "e";
  // append_time(log, 1);
  subspace = make_irrep_subspace(KV_matrix, begin, end, allow_complex);

  // log.indent() << "f";
  // append_time(log, 2);
}

/// Check if Irrep is identity
///
/// - First character is 1+0i, sum of characters == characters.size()
bool PossibleIrrep::is_identity() const {
  std::complex<double> complex_one{1., 0.};
  std::complex<double> first = characters(0);
  std::complex<double> complex_size{double(characters.size()), 0.};
  std::complex<double> sum = characters.sum();

  return almost_equal(first, complex_one, TOL) &&
         almost_equal(sum, complex_size, TOL);
}

/// Check if Irrep is gerade
///
/// - "gerade": inversion results in no sign change
/// - "ungerade": inversion results in sign change
bool PossibleIrrep::is_gerade() const {
  std::complex<double> first = characters(0);
  std::complex<double> last = characters(characters.size() - 1);
  return almost_equal(first, last, TOL);
}

bool PossibleIrrep::operator<(PossibleIrrep const &other) const {
  // Identity comes first
  bool this_is_identity = this->is_identity();
  bool other_is_identity = other.is_identity();
  if (this_is_identity != other_is_identity) {
    return this_is_identity;
  }

  // Low-dimensional irreps come before higher dimensional
  if (!almost_equal(this->characters(0), other.characters(0))) {
    return this->characters(0).real() < other.characters(0).real();
  }

  // 'gerade' irreps come before 'ungerade' irreps
  // This check may need to be improved to know whether inversion is actually
  // present
  bool this_is_gerade = this->is_gerade();
  bool other_is_gerade = other.is_gerade();
  if (this_is_gerade != other_is_gerade) {
    return this_is_gerade;
  }

  // Finally, compare lexicographically (real first, then imag)
  for (Index i = 0; i < this->characters.size(); ++i) {
    if (!almost_equal(this->characters(i).real(), other.characters(i).real()))
      return this->characters(i).real() > other.characters(i).real();
  }
  for (Index i = 0; i < this->characters.size(); ++i) {
    if (!almost_equal(this->characters(i).imag(), other.characters(i).imag()))
      return this->characters(i).imag() > other.characters(i).imag();
  }

  // Now, any possible irrep that are still equal, we break the tie by
  // comparing subspace vectors
  if (this->subspace.cols() != other.subspace.cols()) {
    return this->subspace.cols() < other.subspace.cols();
  }
  for (Index col = 0; col < this->subspace.cols(); ++col) {
    for (Index i = 0; i < this->subspace.size(); ++i) {
      if (!almost_equal(this->subspace(i, col).real(),
                        other.subspace(i, col).real()))
        return this->subspace(i, col).real() > other.subspace(i, col).real();
    }
    for (Index i = 0; i < this->subspace.size(); ++i) {
      if (!almost_equal(this->subspace(i, col).imag(),
                        other.subspace(i, col).imag()))
        return this->subspace(i, col).imag() > other.subspace(i, col).imag();
    }
  }

  throw std::runtime_error("Error comparing PossibleIrrep, tied");

  return false;
}

/// Given kernel, K, and commuter matrix, M, perform eigenvalue decomposition
///     K.adjoint() * M * K = V * D * V.inverse()
/// and construct matrix representation that acts on vectors in the K*V basis,
/// which will be block diagonalized and sorted by eigenvalue. Each block
/// corresponds to a possible irrep, which can be checked by characters value.
std::vector<PossibleIrrep> make_possible_irreps(
    Eigen::MatrixXcd const &commuter, Eigen::MatrixXcd const &kernel,
    MatrixRep const &rep, std::vector<Index> const &head_group_vec,
    bool allow_complex, std::optional<Log> log) {
  // magnify the range of eigenvalues to be (I think) independent of
  // matrix dimension by multiplying to dim^{3/2}
  //
  // solve for eigenvalues and eigenvectors of:
  //    dim^(3/2) * kernel.adjoint() * M * kernel

  //

  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Begin eigenvalue decomposition";
    append_time(*log, 1);
  }

  double dim = kernel.rows();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> esolve;
  double scale = dim * sqrt(dim);
  esolve.compute(scale * kernel.adjoint() * commuter * kernel);
  Eigen::MatrixXd eigenvalues = esolve.eigenvalues();
  Eigen::MatrixXcd KV_matrix = kernel * esolve.eigenvectors();
  Eigen::MatrixXcd KV_adj = KV_matrix.adjoint();

  // Columns of KV_matrix are orthonormal eigenvectors of commuter in terms of
  // natural basis (they were calculated in terms of kernel as basis)

  // When the matrix representation is transformed to operate on coordinates
  // with KV_matrix basis, it becomes block diagonalized
  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Begin block diagonalization";
    append_time(*log, 1);
  }

  // This vector holds the transformed representation matrices for each
  // element indicated by `head_group_vec`, sequentially
  std::vector<Eigen::MatrixXcd> transformed_rep;
  transformed_rep.resize(head_group_vec.size());

  auto worker = [&](Index start, Index end, Index thread_id) {
    Eigen::MatrixXcd temp;
    for (Index idx = start; idx < end; ++idx) {
      Index element_index = head_group_vec[idx];
      temp.noalias() = rep[element_index] * KV_matrix;
      transformed_rep[idx].noalias() = KV_adj * temp;
    }
  };

  threaded_run(head_group_vec.size(), worker);

  // make possible irreps:
  // - The possible irrep corresponds to a range eigenvectors with equal
  //   eigenvalues, could be irrep or could be reducible with degenerate
  //   eigenvalues
  // - When the possible irrep for a range of equal eigenvectors is
  // constructed,
  //   its characters vector is constructed, and if the squared norm of the
  //   characters vectors equals the head group size, then the corresponding
  //   columns of the KV_matrix are an irrep subspace
  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Begin irrep identification... ";
    append_time(*log, 1);
  }

  std::vector<PossibleIrrep> possible_irreps;
  Index begin = 0;
  do {
    double max_eigenvalue = eigenvalues.cwiseAbs().maxCoeff();
    Eigen::VectorXd normalized_eigenvalues = eigenvalues;
    if (max_eigenvalue > 1.0) {
      normalized_eigenvalues /= max_eigenvalue;
    }

    Index end = find_end_of_equal_eigenvalues(begin, normalized_eigenvalues);

    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "- Checking cols [" << begin << ", " << end
                    << ") with eigenvalue = " << std::setprecision(16)
                    << eigenvalues(begin)
                    << " (normalized = " << normalized_eigenvalues(begin)
                    << ")";
      append_time(*log, 1);
    }

    // Should be true by construction for this method
    bool is_block_diagonal = true;

    possible_irreps.emplace_back(eigenvalues, KV_matrix, transformed_rep,
                                 is_block_diagonal, head_group_vec.size(),
                                 allow_complex, begin, end);

    if (log.has_value() && log->verbosity() >= Log::verbose) {
      log->indent() << "  - Is irrep = " << std::boolalpha
                    << possible_irreps.back().is_irrep;
      append_time(*log, 1);
      auto const &p = possible_irreps.back();

      bool is_complex = !almost_zero(p.subspace.adjoint().imag());

      // double max_imag_component = 0.0;
      // if (p.subspace.size() > 0) {
      //   max_imag_component =
      //   p.subspace.adjoint().imag().cwiseAbs().maxCoeff();
      // }

      log->indent() << "    - characters_squared_norm: " << std::setprecision(2)
                    << p.characters_squared_norm << std::endl;
      log->indent() << "    - is_block_diagonal: " << std::boolalpha
                    << p.is_block_diagonal << std::endl;
      log->indent() << "    - complex: " << std::boolalpha << is_complex
                    << std::endl;
      // log->indent() << "    - max imag component of subspace: "
      //               << max_imag_component << std::endl;
      // std::cout << "subspace:\n" << p.subspace << std::endl;
    }

    begin = end;
  } while (begin != eigenvalues.size());

  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Irrep identification: DONE";
    append_time(*log, 2);
  }

  return possible_irreps;
}

/// Make a vector of IrrepInfo from PossibleIrreps
std::vector<IrrepInfo> make_irrep_info(std::set<PossibleIrrep> const &irreps) {
  std::vector<IrrepInfo> irrep_info;
  for (PossibleIrrep const &irrep : irreps) {
    irrep_info.emplace_back(irrep.subspace.adjoint(), irrep.characters);

    if (irrep.subspace.cols() == 2 * (irrep.end - irrep.begin))
      irrep_info.back().pseudo_irrep = true;
    else
      irrep_info.back().pseudo_irrep = false;
  }

  // set sequential indices to differentiate irreps
  // with identical character vectors
  if (irrep_info.size() < 2) {
    return irrep_info;
  }
  Index irrep_index = 0;
  for (Index i = 0; i < irrep_info.size() - 1; ++i) {
    irrep_info[i].index = irrep_index;
    if (almost_equal(irrep_info[i + 1].characters, irrep_info[i].characters,
                     TOL)) {
      irrep_index++;
    } else {
      irrep_index = 0;
    }
  }

  return irrep_info;
}

/// \brief Transforms IrrepInfo constructed for a subspace to be IrrepInfo
/// appropriate for the full space (full space dimension == subspace.rows())
///
/// \param irrep IrrepInfo constructed for a subspace (irrep.trans_mat shape
/// is
///     (subspace.cols() x subspace.rows())
/// \param subspace The subspace that subspace_irrep was constructed for
///
/// \result IrrepInfo constructed for the full space (result.trans_mat shape
/// is
///     (subspace.cols() x subspace.cols())
///
IrrepInfo subspace_to_full_space(IrrepInfo const &subspace_irrep,
                                 Eigen::MatrixXd const &subspace) {
  IrrepInfo result(subspace_irrep);

  result.trans_mat = subspace_irrep.trans_mat *
                     subspace.adjoint().template cast<std::complex<double>>();

  result.irrep_dim = result.trans_mat.rows();
  result.vector_dim = result.trans_mat.cols();

  result.directions.clear();
  for (const auto &direction_orbit : subspace_irrep.directions) {
    std::vector<Eigen::VectorXd> new_orbit;
    new_orbit.reserve(direction_orbit.size());
    for (const auto &directions : direction_orbit) {
      new_orbit.push_back(subspace * directions);
    }
    result.directions.push_back(std::move(new_orbit));
  }
  return result;
}

/// Check if a representation is irreducible
///
/// A representation is irreducible if the squared norm of the characters
/// equals the group size
bool is_irrep(MatrixRep const &rep, GroupIndices const &head_group) {
  double characters_squared_norm = 0;
  for (Index element_index : head_group) {
    double character = rep[element_index].trace();
    characters_squared_norm += character * character;
  }
  return almost_equal(characters_squared_norm, double(head_group.size()), TOL);
}

/// IrrepDecomposition proceeds by constructing "commuters", M_k, which
/// commute (M_k * R(r) = R(r) * M_k) with all of the matrix representations,
/// R(r), of the group. The commuters are constructed to reveal irreducible
/// vector spaces (via application of a Reynolds operator), and be orthonormal
/// to existing commuters (via Gram-Shmidt). The commuters are found via a
/// process which constructs a candidate commuter which is either the Zero
/// matrix, and then skipped, or else it is a useful non-zero commuter which
/// will block diagonalize :
///
///     M_candidate(i,j,phase) = sum_r R(r) * M_init * R(r).transpose()
///     M_init = phase * K.col(i) * K.col(j).adjoint() +
///              std::conj(phase) * K.col(j) * K.col(i).adjoint()
///
/// where:
/// - K: kernel matrix, the null space of the already found irreducible vector
/// spaces. The kernel matrix is initialized as a full rank matrix and over
/// the course of the IrrepDecomposition the kernel shrinks as the irreducible
/// vector spaces are found.
/// - candidate commuting matrices, M_candidate, are built from the outer
/// product of two columns, i, and j, of the kernel matrix, and a complex
/// phase parameter (1 or i)
/// - R(r): is the matrix representation for element r of the head_group
/// - M_k: previously found commuting matrices
///
/// Once a new non-zero commuter is found, possible irreducible subspaces are
/// found and checked. A possible irreducible subspace is each
/// K*V_equal_eigenvalue_set[i], where V_equal_eigenvalue_set[i] is the vector
/// space corresponding to eigenvectors of (K.adjoint() * M_new * K) with
/// equal eigenvalues.
///
/// Eigenvalue decomposition:
///
///     K.adjoint() * M_new * K = V * D * V.inverse()
///
/// - D: diagonal matrix of sorted eigenvalues (size K.cols() x K.cols())
/// - V: eigenvector column matrix (size K.cols() x K.cols())
///
/// For each set of equal eigenvalues, a PossibleIrrep is constructed that
/// stores:
/// - `begin`: the column of the first of the set of equal eigenvalues
/// - `irrep_dim`: the subspace dimension / the number of equal eigenvalues
/// - `subspace`: irreducible subspace, (dim x irrep_dim (?) matrix):
///
///   First step, find vector space:
///   - If allow_complex: subspace = X
///   - If !allow_complex: subspace =
///         [sqrt(2.0) * X.real(), sqrt(2.0) *  X.imag()],
///     where X = (K * V).block(0, begin,
///                             K.rows(), irrep_dim)
///   Second step, orthogonalize via QR decomposition
///
/// - `characters`: Vector of complex characters. The character of a matrix
///   representation is the trace of the representation.
/// - `characters_squared_norm`: For an irreducible representation, the
/// squared
///   norm of the characters vector is equal the size of the group.
/// - `symmetrizer`: A pair with, symmetrizer.first being a MatrixXcd, which
///   defines a rotation of the irreducible subspace that aligns its
///   components along high-symmetry directions, and symmetrizer.second being
///   a vector of orbits of high-symmetry directions
///
/// For each PossibleIrrep, check if it extends adapted_subspace. If it
/// does, then symmetrize and save the irrep. For all new irreps, extend
/// adapted_subspace to include the irrep's subspace. Then once all the new
/// irreps are added, recalculate the kernel matrix and begin the loop again,
/// until the adapted_subspace is full rank.

/// Finds irreducible subspaces that comprise an underlying subspace
///
/// This method does not rely on the character table, but instead utilizes a
/// brute-force approach. It is not guaranteed to find all irreps, so the
/// resulting irreps should be checked if they span the entire space
/// represented by `rep`. This can be done by checking if
///     `full_trans_mat(result).adjoint().rows() == rep[i].rows()`).
/// This method does not align the irrep subspace axes along high symmetry
/// directions.
///
/// \param rep Matrix representation of head_group, this defines group action
/// on the underlying vector space
/// \param head_group Group for which the irreps are to be found
/// \param allow_complex If true, irreducible space basis vectors may be
///     complex-valued. If false, complex irreps are combined to form real
///     representations
/// \param log Optional Log object for logging progress
/// \param method Method for constructing commuters. The "deterministic" method
///     iterates deterministically through possible commuters, but does not
///     guarantee that all irreps will be found. The "random" method randomly
///     randomly constructs commuters, alternating between commuters
///     constructed from real-valued and complex-valued seeds, up to a maximum
///     of 10 attempts.
/// \param start_with_real_seed If method == CommuterMethod::random, whether to
///     start with a commuter constructed from a real-valued seed, or a
///     complex-valued seed.
///
/// \result vector of IrrepInfo objects. Irreps are ordered by dimension, with
///     identity first (if present).  Repeated irreps (with equal character
///     vectors) are sequential, and are distinguished by IrrepInfo::index.
///
std::vector<IrrepInfo> irrep_decomposition(
    MatrixRep const &rep, GroupIndices const &head_group, bool allow_complex,
    std::optional<Log> log, CommuterMethod method, bool start_with_real_seed) {
  if (log.has_value()) {
    log->begin<Log::standard>("Find irreps");
    log->increase_indent();
    log->indent() << std::endl;

    log->indent() << "Using " << max_threads() << " threads" << std::endl
                  << std::endl;
    log->indent() << "Number of group elements = " << rep.size() << std::endl;
  }

  if (!rep.size()) {
    if (log.has_value()) {
      log->indent() << std::endl;
      log->indent() << "No irreps to find." << std::endl << std::endl;
      log->decrease_indent();
    }
    return std::vector<IrrepInfo>();
  }

  std::vector<Index> head_group_vec;
  head_group_vec.reserve(head_group.size());
  for (Index idx : head_group) head_group_vec.push_back(idx);

  if (log.has_value()) {
    log->indent() << "Vector space dimension = " << rep[0].rows() << std::endl
                  << std::endl;
  }

  int dim = rep[0].rows();

  // This method iteratively finds irreducible spaces, which are used to
  // extend the "adapted_subspace" (combined space of found irreducible
  // spaces). The "adapted_subspace" is not aligned along high symmetry
  // directions by this function. When the "adapted_subspace" is of dimenions
  // equal to `dim`, all irreps have been found.

  // start with all kernel, end with all adapted_subspace
  Eigen::MatrixXcd kernel = complex_I(dim, dim);
  Eigen::MatrixXcd adapted_subspace{dim, 0};

  // In this set, as they are discovered we will save PossibleIrrep that:
  // - i) actually are irreducible,
  // - and ii) have distinct subspaces (detected when their subspace extends
  //   the adapted_subspace space) BP: not necessary?
  std::set<PossibleIrrep> irreps;

  if (method == CommuterMethod::deterministic) {
    // count over possible commuter matrices for this kernel:
    // - kernel column pairs (i,j), j>=i & phase = [1, i]; skips i==j if
    // phase==i
    CommuterParamsCounter commuter_params;
    commuter_params.reset(kernel);

    do {  // while adapated_subspace.cols() != dim

      if (!commuter_params.valid()) {
        // The commuter construction method does not currently guarantee that
        // all irreps will be revealed. The caller may have a way to handle
        // this and so this does not throw an exception.

        if (log.has_value()) {
          log->indent() << std::endl;
          log->indent() << "Break: All commuters attempted" << std::endl;
        }

        break;
      }
      if (log.has_value() && log->verbosity() >= Log::verbose) {
        log->indent() << "Make commuter... ";
        append_time(*log, 1);
      }

      // make next commuter, M, and check if not zero
      Eigen::MatrixXcd commuter =
          make_commuter(commuter_params, rep, head_group, kernel);

      if (almost_equal(frobenius_product(commuter).real(), 0., TOL)) {
        if (log.has_value() && log->verbosity() >= Log::verbose) {
          log->indent() << "Frobenius product is zero. Skipping... ";
          append_time(*log, 1);
        }
        commuter_params.increment();
        continue;
      }

      // make possible irreps:
      //
      // Given kernel, K, and commuter matrix, M, perform eigenvalue
      // decomposition
      //     K.adjoint() * M * K = V * D * V.inverse()
      // and construct matrix representation that acts on vectors in the K*V
      // basis, which will be block diagonalized and sorted by eigenvalue. Each
      // block corresponds to a possible irrep, which can be checked by its
      // characters. The columns in K*V corresponding to an irrep are the irrep
      // subspace.
      if (log.has_value() && log->verbosity() >= Log::verbose) {
        log->indent() << "Make possible irreps...";
        append_time(*log, 1);
      }
      std::vector<PossibleIrrep> possible_irreps = make_possible_irreps(
          commuter, kernel, rep, head_group_vec, allow_complex, log);

      // save any possible irrep that:
      // - i) is an irrep,
      // - and ii) extends the adapted_subspace space (BP: not necessary?)
      bool any_new_irreps = false;
      for (auto const &possible_irrep : possible_irreps) {
        if (possible_irrep.is_irrep &&
            is_extended_by(adapted_subspace, possible_irrep.subspace)) {
          irreps.insert(possible_irrep);
          adapted_subspace = extend(adapted_subspace, possible_irrep.subspace);
          any_new_irreps = true;

          if (log.has_value()) {
            log->indent() << "Found irrep of dim " << possible_irrep.irrep_dim
                          << " (" << dim - adapted_subspace.cols() << " / "
                          << dim << " dim remaining)";
            append_time(*log, 1);
          }
        }
      }

      // If any new irreps were found, break to return
      // Empirically, it seems more efficient to break and continue
      // in the smaller remaining subspace with recalculated matrix reps.
      if (any_new_irreps && adapted_subspace.cols() != dim) {
        kernel = make_kernel(adapted_subspace);
        commuter_params.reset(kernel);
        if (kernel.cols() + adapted_subspace.cols() !=
            adapted_subspace.rows()) {
          throw std::runtime_error(
              "Unknown error finding irreps: dimension mismatch");
        }
        // Restart after finding any irreps
        if (log.has_value()) {
          log->indent() << std::endl;
          log->indent() << "Break: irreps found" << std::endl;
        }

        break;

      } else {
        commuter_params.increment();
      }
    } while (adapted_subspace.cols() != dim);

  } else {
    // CommuterMethod::random
    // Use random seed matrices projected via Reynolds operator.
    // Alternates between real symmetric and complex Hermitian seeds:
    // - Real symmetric seeds find real irreps with clean real eigenvectors
    // - Complex Hermitian seeds can distinguish complex conjugate irrep
    //   pairs that share eigenvalues under a real symmetric commuter
    std::random_device rd;
    std::mt19937 gen(rd());
    Index max_random_attempts = 10;

    bool any_new_irreps = false;
    Index attempt = 0;

    while (attempt < max_random_attempts) {
      // Alternate real seed and complex seed
      Index offset = start_with_real_seed ? 0 : 1;
      bool use_complex_seed = ((offset + attempt) % 2 == 1);

      if (log.has_value() && log->verbosity() >= Log::verbose) {
        log->indent() << "Make random commuter (attempt " << attempt + 1
                      << " / " << max_random_attempts << ", "
                      << (use_complex_seed ? "complex" : "real")
                      << " seed)... ";
        append_time(*log, 1);
      }

      Eigen::MatrixXcd commuter =
          make_random_commuter(rep, head_group, kernel, gen, use_complex_seed);

      if (almost_equal(frobenius_product(commuter).real(), 0., TOL)) {
        if (log.has_value() && log->verbosity() >= Log::verbose) {
          log->indent() << "Frobenius product is zero. Skipping... ";
          append_time(*log, 1);
        }
        ++attempt;
        continue;
      }

      if (log.has_value() && log->verbosity() >= Log::verbose) {
        log->indent() << "Make possible irreps...";
        append_time(*log, 1);
      }
      std::vector<PossibleIrrep> possible_irreps = make_possible_irreps(
          commuter, kernel, rep, head_group_vec, allow_complex, log);

      for (auto const &possible_irrep : possible_irreps) {
        if (possible_irrep.is_irrep &&
            is_extended_by(adapted_subspace, possible_irrep.subspace)) {
          irreps.insert(possible_irrep);
          adapted_subspace = extend(adapted_subspace, possible_irrep.subspace);
          any_new_irreps = true;

          if (log.has_value()) {
            log->indent() << "Found irrep of dim " << possible_irrep.irrep_dim
                          << " (" << dim - adapted_subspace.cols() << " / "
                          << dim << " dim remaining)";
            append_time(*log, 1);
          }
        }
      }

      if (any_new_irreps) {
        break;
      }
      ++attempt;
    }

    if (log.has_value()) {
      if (!any_new_irreps) {
        log->indent() << std::endl;
        log->indent() << "Break: No new irreps found after "
                      << max_random_attempts << " random commuter attempts"
                      << std::endl;
      } else {
        log->indent() << std::endl;
        log->indent() << "Break: irreps found" << std::endl;
      }
    }
  }

  // Make irrep info (no directions yet, orthogonalized but not aligned along
  // high symmetry directions)
  std::vector<IrrepInfo> irrep_info = make_irrep_info(irreps);

  if (log.has_value()) {
    log->indent() << std::endl;
    log->indent() << "Found " << irrep_info.size() << " irreps." << std::endl;
    log->indent() << "Found irreps for " << adapted_subspace.cols() << " / "
                  << dim << " dimensions." << std::endl;
    log->indent() << "Complete irrep decomposition: "
                  << (adapted_subspace.cols() == dim ? "yes" : "no")
                  << std::endl
                  << std::endl;
    log->decrease_indent();
  }
  return irrep_info;
}

/// Convert irreps generated for a subspace to full space dimension
///
/// \param subspace_irreps Irreducible spaces in the subspace
/// (subspace_irreps[i].trans_mat.rows() == subspace dimension,
/// subspace_irreps[i].trans_mat.cols() == fullspace dimensino) \param
/// subspace Basis for a subspace (subspace.rows() == fullspace
///     dimension, subspace.cols() == subspace dimension)
std::vector<IrrepInfo> make_fullspace_irreps(
    std::vector<IrrepInfo> const &subspace_irreps,
    Eigen::MatrixXd const &subspace) {
  std::vector<IrrepInfo> fullspace_irreps;
  fullspace_irreps.reserve(subspace_irreps.size());
  for (auto const &irrep : subspace_irreps) {
    fullspace_irreps.push_back(subspace_to_full_space(irrep, subspace));
  }
  return fullspace_irreps;
}

/// Find the invariant subspace generated by applying group to subspace
/// using modified Gram-Schmidt (multithreaded)
Eigen::MatrixXd make_invariant_space(MatrixRep const &rep,
                                     GroupIndices const &head_group,
                                     Eigen::MatrixXd const &subspace) {
  std::vector<Index> head_group_vec;
  head_group_vec.reserve(head_group.size());
  for (Index idx : head_group) head_group_vec.push_back(idx);
  return make_invariant_space(rep, head_group_vec, subspace);
}

/// Find the invariant subspace generated by applying group to subspace
/// using modified Gram-Schmidt (multithreaded)
Eigen::MatrixXd make_invariant_space(MatrixRep const &rep,
                                     std::vector<Index> const &head_group_vec,
                                     Eigen::MatrixXd const &subspace) {
  const double tol = TOL;
  Index n = subspace.rows();
  Index k = subspace.cols();

  // If there are no group elements, return empty basis
  if (head_group_vec.size() == 0 || k == 0) {
    return Eigen::MatrixXd(n, 0);
  }

  Index n_threads = max_threads();

  // Prepare per-thread local bases
  std::vector<Eigen::MatrixXd> local_bases;
  local_bases.resize(n_threads);

  auto worker = [&](Index start, Index end, Index t) {
    Eigen::MatrixXd &basis = local_bases[t];
    basis.resize(n, 0);

    for (Index idx = start; idx < end; ++idx) {
      Index element_index = head_group_vec[idx];
      Eigen::MatrixXd transformed = rep[element_index] * subspace;  // (n x k)

      for (Index c = 0; c < k; ++c) {
        Eigen::VectorXd v = transformed.col(c);
        // Modified Gram-Schmidt against thread-local basis
        for (Index j = 0; j < basis.cols(); ++j) {
          double proj = basis.col(j).dot(v);
          v.noalias() -= proj * basis.col(j);
        }

        double norm = v.norm();
        if (norm > tol) {
          v /= norm;
          Index old_cols = basis.cols();
          basis.conservativeResize(n, old_cols + 1);
          basis.col(old_cols) = std::move(v);
        }
      }
    }
  };

  threaded_run(head_group_vec.size(), worker);

  // Combine local bases into a single global basis using Modified
  // Gram-Schmidt
  Eigen::MatrixXd basis(n, 0);
  for (Index t = 0; t < n_threads; ++t) {
    Eigen::MatrixXd const &lb = local_bases[t];
    for (Index col = 0; col < lb.cols(); ++col) {
      Eigen::VectorXd v = lb.col(col);

      // MGS against current global basis
      for (Index j = 0; j < basis.cols(); ++j) {
        double proj = basis.col(j).dot(v);
        v.noalias() -= proj * basis.col(j);
      }

      double norm = v.norm();
      if (norm > tol) {
        v /= norm;
        Index old_cols = basis.cols();
        basis.conservativeResize(n, old_cols + 1);
        basis.col(old_cols) = std::move(v);
      }
    }
  }

  if (basis.cols() == 0) {
    return Eigen::MatrixXd(n, 0);
  }

  // Clean up near-linear dependencies / ensure orthonormal columns
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> colqr(basis);
  colqr.setThreshold(tol);
  Eigen::MatrixXd Q = colqr.householderQ();
  return Q.leftCols(colqr.rank());
}

/// \brief Create the subspace rep from the fullspace rep
///
/// Create `subspace_rep`, a transformed copy of `fullspace_rep` that acts
/// on coordinates with `subspace` columns as a basis. Matrices in
/// `subspace_rep` are shape (subspace.cols() x subspace.cols())
///
/// Notes: This function uses threads to parallelize the construction of the
/// subspace representation matrices.
///
//// \param fullspace_rep Matrix representation for transforming unrolled
/// vectors in the prim basis
/// \param subspace A subspace basis, x_fullspace = subspace * x_subspace.
/// Subspace basis vectors must be orthonormal.
///
/// \return subspace_rep, The matrix representation for transforming vectors
/// in the subspace
MatrixRep make_subspace_rep(MatrixRep const &fullspace_rep,
                            Eigen::MatrixXd const &subspace) {
  // x_f' = M_f * x_f
  // x_f = B * x_s
  // where:
  // - x_f: vector in full space
  // - M_f: full space matrix representation
  // - x_s: vector in subspace
  // - B: subspace basis matrix (fullspace_dim x subspace_dim)

  // therefore, the subspace representation is:
  //     x_s' = M_s * x_s,
  //     M_s = B_pinv * M_f * B

  // this function assumes B is orthonormal, so B_pinv = B.transpose()

  Index const n = static_cast<Index>(fullspace_rep.size());
  MatrixRep subspace_rep;
  if (n == 0) return subspace_rep;

  // Pre-size result to avoid reallocations during parallel writes
  subspace_rep.resize(n);

  // Define the worker lambda that fills subspace_rep for indices [start,end)
  auto worker = [&](Index start, Index end, Index /*thread_id*/) {
    Eigen::MatrixXd temp;
    for (Index i = start; i < end; ++i) {
      temp.noalias() = fullspace_rep[i] * subspace;
      subspace_rep[i].noalias() = subspace.transpose() * temp;
    }
  };

  threaded_run(n, worker);

  return subspace_rep;
}

/// \brief Symmetrize IrrepInfo, by finding high symmetry directions and
/// aligning the irrep subspace basis with those directions
std::vector<IrrepInfo> symmetrize_irreps(
    MatrixRep const &subspace_rep, GroupIndices const &head_group,
    std::vector<IrrepInfo> const &irreps,
    GroupIndicesOrbitSet const &subgroup_orbits, std::optional<Log> log) {
  std::vector<IrrepInfo> symmetrized_irreps;
  double vec_compare_tol = TOL;

  Index i_irrep = 1;
  for (const auto &irrep : irreps) {
    if (log.has_value() && log->print()) {
      if (log->verbosity() >= Log::verbose) {
        std::stringstream ss;
        ss << "Symmetrize irrep " << i_irrep << " / " << irreps.size();
        log->begin<Log::standard>(ss.str());
        log->indent() << std::endl;
        log->indent() << "Irrep dim = " << irrep.irrep_dim << std::endl;
      } else {
        log->indent() << "Irrep " << i_irrep << " / " << irreps.size()
                      << ":  dim = " << irrep.irrep_dim;
      }
    }

    Eigen::MatrixXcd irrep_subspace = irrep.trans_mat.adjoint();

    multivector<Eigen::VectorXcd>::X<2> irrep_special_directions =
        make_irrep_special_directions(subspace_rep, head_group, irrep_subspace,
                                      vec_compare_tol, subgroup_orbits, log);

    Eigen::MatrixXcd symmetrizer_matrix = make_irrep_symmetrizer_matrix(
        irrep_special_directions, irrep_subspace, vec_compare_tol, log);

    IrrepInfo symmetrized_irrep{irrep};
    symmetrized_irrep.trans_mat =
        (irrep_subspace * symmetrizer_matrix).adjoint();
    symmetrized_irrep.directions = to_real(irrep_special_directions);
    symmetrized_irreps.push_back(symmetrized_irrep);

    if (log.has_value() && log->print()) {
      if (log->verbosity() < Log::verbose) {
        log->indent() << "  orbits of special directions = "
                      << irrep_special_directions.size() << " ";
        append_time(*log, 1);
      }
    }

    i_irrep++;
  }

  if (log.has_value() && log->print()) {
    if (log->verbosity() >= Log::verbose) {
      log->indent() << std::endl;
    }
  }
  return symmetrized_irreps;
}

}  // namespace IrrepDecompositionImpl

}  // namespace irreps

}  // namespace CASM
