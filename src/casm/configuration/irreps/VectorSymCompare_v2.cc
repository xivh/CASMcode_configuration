#include "casm/configuration/irreps/VectorSymCompare_v2.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

namespace irreps {

VectorInvariants::VectorInvariants(Eigen::VectorXcd const &vector)
    : m_cols(vector.cols()), m_norm(vector.norm()) {}

double VectorInvariants::cols() const { return m_cols; }

double VectorInvariants::norm() const { return m_norm; }

}  // namespace irreps

bool almost_equal(irreps::VectorInvariants const &A_invariants,
                  irreps::VectorInvariants const &B_invariants, double tol) {
  return A_invariants.cols() == B_invariants.cols() &&
         almost_equal(A_invariants.norm(), B_invariants.norm(), tol);
}

bool compare(irreps::VectorInvariants const &A_invariants,
             irreps::VectorInvariants const &B_invariants, double tol) {
  if (A_invariants.cols() == B_invariants.cols()) {
    return CASM::compare(A_invariants.norm(), B_invariants.norm(), tol);
  }
  return A_invariants.cols() < B_invariants.cols();
}

namespace irreps {

VectorSymCompare::VectorSymCompare(MatrixRep const &matrix_rep, double tol)
    : m_matrix_rep(matrix_rep), m_tol(tol) {}

// lexicographical comparison (reversed, uses vector_B < vector_A)
bool VectorSymCompare::compare(Eigen::VectorXcd const &vector_A,
                               Eigen::VectorXcd const &vector_B) const {
  return colmajor_lex_compare(vector_B, vector_A, m_tol);
}

// return VectorInvariants
VectorInvariants VectorSymCompare::make_invariants(
    Eigen::VectorXcd const &vector) const {
  return VectorInvariants{vector};
}

// apply matrix rep to vector
Eigen::VectorXcd VectorSymCompare::copy_apply(Index const &op_index,
                                              Eigen::VectorXcd vector) const {
  return m_matrix_rep[op_index] * vector;
}

// no change needed to prepare for comparison
Eigen::VectorXcd VectorSymCompare::prepare(Eigen::VectorXcd vector) const {
  return vector;
}

// compare orbits by first comparing invariants, then orbit prototypes
bool VectorSymCompare::inter_orbit_compare(
    Eigen::VectorXcd const &A_prototype, VectorInvariants const &A_invariants,
    Eigen::VectorXcd const &B_prototype,
    VectorInvariants const &B_invariants) const {
  // first compare invariants
  if (CASM::compare(A_invariants, B_invariants, m_tol)) {
    return true;
  }
  if (CASM::compare(B_invariants, A_invariants, m_tol)) {
    return false;
  }

  // next compare A and B
  return this->compare(A_prototype, B_prototype);
}

namespace {

template <typename Derived>
typename Derived::PlainObject standardize_column_vector_signs_impl(
    Eigen::MatrixBase<Derived> const &vector_space, double tol) {
  typename Derived::PlainObject result(vector_space);
  for (CASM::Index row = 0; row < result.rows(); ++row) {
    for (CASM::Index col = 0; col < result.cols(); ++col) {
      if (!CASM::almost_zero(result(row, col), tol)) {
        result.col(col) *= std::abs(result(row, col)) / result(row, col);
        break;
      }
    }
  }
  return result;
}

template <typename Derived>
typename Derived::PlainObject vector_space_prepare_impl(
    Eigen::MatrixBase<Derived> const &obj, double _tol) {
  typename Derived::PlainObject result =
      typename Derived::PlainObject(CASM::reduced_column_echelon(obj, _tol)
                                        .householderQr()
                                        .householderQ())
          .leftCols(obj.cols());
  return standardize_column_vector_signs_impl(result, _tol);
}
}  // namespace

/// Standardize basis vector signs
Eigen::MatrixXcd standardize_column_vector_signs(
    Eigen::MatrixXcd const &vector_space, double tol) {
  return standardize_column_vector_signs_impl(vector_space, tol);
}

/// Standardize basis vector signs
Eigen::MatrixXd standardize_column_vector_signs(
    Eigen::MatrixXd const &vector_space, double tol) {
  return standardize_column_vector_signs_impl(vector_space, tol);
}

/// Vector space preparation for comparison
///
/// - Attempts to find a sparse set of spanning vectors, and sort them so that
/// subspace matrix is nearly upper triangular (if possible)
/// - Also enures that first nonzero element of each row (if there is one) is
/// positive
Eigen::MatrixXcd vector_space_prepare(Eigen::MatrixXcd const &vector_space,
                                      double tol) {
  return vector_space_prepare_impl(vector_space, tol);
}

/// Vector space preparation for comparison
///
/// - Attempts to find a sparse set of spanning vectors, and sort them so that
/// subspace matrix is nearly upper triangular (if possible)
/// - Also enures that first nonzero element of each row (if there is one) is
/// positive
Eigen::MatrixXd vector_space_prepare(Eigen::MatrixXd const &vector_space,
                                     double tol) {
  return vector_space_prepare_impl(vector_space, tol);
}

}  // namespace irreps

}  // namespace CASM
