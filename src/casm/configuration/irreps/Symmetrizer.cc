#include "casm/configuration/irreps/Symmetrizer.hh"

#include "casm/casm_io/container/stream_io.hh"
#include "casm/configuration/irreps/SimpleOrbit_impl.hh"
#include "casm/configuration/irreps/VectorSymCompare_v2.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

// debug
#include <iostream>

namespace CASM {

namespace irreps {

/// Find high-symmetry directions in a irreducible space
///
/// \param rep Matrix representation of head_group, this defines group action
/// on the underlying vector space
/// \param head_group Group for which the irreps are to be found
/// \param subspace A column vector matrix representing a basis of the
///     irreducible space in which high symmetry directions will be found. The
///     number of rows must equal `rep.dim()`, the number of columns is equal to
///     the dimension of the irreducible space.
/// \param vec_compare_tol Tolerance for elementwise floating-point comparisons
///     of vectors
/// \param all_subgroups Denotes whether all subgroups of head_group should be
///     used for symmetry analysis (if true), or only cyclic subgroups (if
///     false). Cyclic subgroups are those found by taking a group element and
///     multiplying it by itself until a group is generated.
///
/// \result Set of directions in the vector space on which 'rep' is defined,
/// such that each direction is invariant to a unique subgroup of 'head_group'
/// (i.e., no other direction in the space, except the negative of that
/// direction, is invariant to that subgroup). The value `result[i]` is an
/// orbit of symmetrically equivalent directions, and the value `result[i][j]`
/// is an individual direction (Eigen::VectorXcd). Direction vectors are
/// normalized to unit length. The total set of all directions is guaranteed to
/// span the space.
///
/// \throws if `rep` is not an irreducible representation
///
multivector<Eigen::VectorXcd>::X<2> make_irrep_special_directions(
    MatrixRep const &rep, GroupIndices const &head_group,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol,
    std::function<GroupIndicesOrbitSet()> make_subgroups_f,
    std::optional<Log> log) {
  if (log.has_value()) {
    log->indent() << "Get subgroup indices..." << std::endl;
  }

  GroupIndicesOrbitSet sgroups = make_subgroups_f();

  if (log.has_value()) {
    log->indent() << "Get subgroup indices: DONE" << std::endl << std::endl;
    log->indent() << "Number of subgroups = " << sgroups.size() << std::endl
                  << std::endl;
    log->indent() << "Applying Reynolds operator to find special directions... "
                  << std::endl
                  << std::endl;
    log->indent() << "Special directions: ";
  }

  std::vector<Eigen::VectorXcd> tdirs;
  Eigen::MatrixXd R;
  Index dim = rep[0].rows();

  VectorSymCompare sym_compare{rep, vec_compare_tol};
  std::set<SimpleOrbit<VectorSymCompare>> orbit_result;

  // Loop over small (i.e., cyclic) subgroups and hope that each special
  // direction is invariant to at least one small subgroup
  for (auto const &orbit : sgroups) {
    // Reynolds for small subgroup *(orbit.begin()) in irrep_subspace i
    R.setZero(dim, dim);

    for (Index element_index : *(orbit.begin())) {
      R += rep[element_index];
    }

    if ((R * irrep_subspace).norm() < TOL) continue;

    // Find spanning vectors of column space of R*irrep_space, which is
    // projection of irrep_space into its invariant component
    auto QR = (R * irrep_subspace).colPivHouseholderQr();
    QR.setThreshold(TOL);

    // If only one spanning vector, it is special direction
    if (QR.rank() > 1) {
      continue;
    }
    Eigen::MatrixXcd Q = QR.matrixQ();

    // Convert from irrep_subspace back to total space and push_back
    tdirs.push_back(Q.col(0));
    auto result_1 = orbit_result.emplace(Q.col(0), head_group.begin(),
                                         head_group.end(), sym_compare);
    if (result_1.second == true) {
      if (log.has_value() && log->print()) {
        log->ostream() << "*";
      }
    }

    tdirs.push_back(-Q.col(0));
    auto result_2 = orbit_result.emplace(-Q.col(0), head_group.begin(),
                                         head_group.end(), sym_compare);
    if (result_2.second == true) {
      if (log.has_value() && log->print()) {
        log->ostream() << "*";
      }
    }
  }

  if (log.has_value() && log->print()) {
    log->ostream() << std::endl << std::endl;
  }

  // t_result may contain duplicates, or elements that are equivalent by
  // symmetry. To discern more info, we need to exclude duplicates and find
  // the orbit of the directions. this should also
  // reveal the invariant subgroups.

  // VectorSymCompare sym_compare{rep, vec_compare_tol};
  // std::set<SimpleOrbit<VectorSymCompare>> orbit_result;
  // for (Eigen::VectorXcd const &direction : tdirs) {
  //   orbit_result.emplace(direction, head_group.begin(), head_group.end(),
  //                        sym_compare);
  // }
  multivector<Eigen::VectorXcd>::X<2> result;
  for (auto const &orbit : orbit_result) {
    result.emplace_back(orbit.begin(), orbit.end());
  }

  if (log.has_value() && log->print()) {
    log->indent() << "Found " << result.size()
                  << " orbits of special directions." << std::endl
                  << std::endl;
  }

  return result;
}

/// Make an irreducible space symmetrizer matrix using special directions
Eigen::MatrixXcd make_irrep_symmetrizer_matrix(
    multivector<Eigen::VectorXcd>::X<2> const &irrep_special_directions,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol,
    std::optional<Log> log) {
  // Four strategies, in order of desparation
  // 1) find a spanning set of orthogonal axes within a single orbit of special
  // directions 2) find a spanning set of orthogonal axes within the total set
  // of special directions 3) perform qr decomposition on lowest-multiplicity
  // orbit to find spanning set of axes 4) find column-echelon form of
  // irrep_subspace matrix to get a sparse/pretty set of axes
  Eigen::MatrixXcd result;
  Index dim = irrep_subspace.cols();
  Index min_mult = 10000;
  bool orb_orthog = false;
  bool tot_orthog = false;

  Index i_strategy = -1;
  Index i_orb_orthog = -1;
  Index i_min_mult_orbit = -1;
  std::vector<Index> i_tot_orthog;
  Eigen::MatrixXcd axes, orb_axes, tot_axes;
  Index tot_col(0);
  tot_axes.setZero(irrep_subspace.rows(), dim);

  // std::cout << "BEGIN MAKE SYMMETRIZED AXES" << std::endl;
  if (irrep_special_directions.size() && log.has_value()) {
    log->indent() << "Find symmetrized axes using special directions..."
                  << std::endl;
  }
  Index i_orbit = 0;
  for (auto const &orbit : irrep_special_directions) {
    // Strategy 1
    if ((orb_orthog && orbit.size() < min_mult) || !orb_orthog) {
      orb_axes.setZero(irrep_subspace.rows(), dim);
      Index col = 0;
      for (auto const &el : orbit) {
        if (almost_zero((el.adjoint() * orb_axes).eval(), vec_compare_tol)) {
          if (col < orb_axes.cols())
            orb_axes.col(col++) = el;
          else {
            std::stringstream errstr;
            errstr
                << "Error in irrep_symmtrizer_from_directions(). Constructing "
                   "coordinate axes from special directions of space spanned "
                   "by row vectors:\n "
                << irrep_subspace.transpose()
                << "\nAxes collected thus far are the row vectors:\n"
                << orb_axes.transpose()
                << "\nBut an additional orthogonal row vector has been found:\n"
                << el.transpose()
                << "\nindicating that irrep_subspace matrix is malformed.";
            throw std::runtime_error(errstr.str());
          }
        }
      }
      if (col == dim) {
        i_strategy = 1;
        i_orb_orthog = i_orbit;
        // std::cout << "PLAN A-- col: " << col << "; dim: " << dim << ";
        // min_mult: " << min_mult << "; orthog: " << orb_orthog << ";\naxes:
        // \n"
        // << axes << "\norb_axes: \n" << orb_axes << "\n\n";
        orb_orthog = true;
        min_mult = orbit.size();
        axes = orb_axes;
      }
    }

    // Greedy(ish) implementation of strategy 2 -- may not find a solution, even
    // if it exists
    if (!orb_orthog && !tot_orthog) {
      for (auto const &el : orbit) {
        if (almost_zero((el.adjoint() * tot_axes).eval(), vec_compare_tol)) {
          if (tot_col < tot_axes.cols()) {
            // std::cout << "Strategy 2: added axis from orbit=" << i_orbit
            //           << std::endl;
            tot_axes.col(tot_col++) = el;
            i_tot_orthog.push_back(i_orbit);
          } else {
            std::stringstream errstr;
            errstr
                << "Error in irrep_symmtrizer_from_directions(). Constructing "
                   "coordinate axes from special directions of space spanned "
                   "by row vectors:\n "
                << irrep_subspace.transpose()
                << "\nAxes collected thus far are the row vectors:\n"
                << tot_axes.transpose()
                << "\nBut an additional orthogonal row vector has been found:\n"
                << el.transpose()
                << "\nindicating that irrep_subspace matrix is malformed.";
            throw std::runtime_error(errstr.str());
          }
        }
      }
      if (tot_col == dim) {
        // std::cout << "Success: Strategy 2" << std::endl;
        // std::cout << "PLAN B-- col: " << tot_col << "; dim: " << dim << ";
        // min_mult: " << min_mult << "; orthog: " << tot_orthog << ";\naxes:
        // \n"
        // << axes << "\ntot_axes: \n" << tot_axes << "\n\n";
        i_strategy = 2;
        tot_orthog = true;
        axes = tot_axes;
      }
    }

    // Strategy 3
    if (!orb_orthog && !tot_orthog && orbit.size() < min_mult) {
      orb_axes.setZero(irrep_subspace.rows(), orbit.size());
      for (Index col = 0; col < orbit.size(); ++col) {
        orb_axes.col(col) = orbit[col];
      }
      // std::cout << "Strategy 3: orbit=" << i_orbit << " mult=" <<
      // orbit.size()
      //           << std::endl;
      // std::cout << "PLAN C--  dim: " << dim << "; min_mult: " << min_mult <<
      // "; orthog: " << orb_orthog << ";\naxes: \n" << axes;
      min_mult = orbit.size();
      axes = Eigen::MatrixXcd(orb_axes.colPivHouseholderQr().matrixQ())
                 .leftCols(dim);
      i_min_mult_orbit = i_orbit;
      i_strategy = 3;
      // std::cout << "\norb_axes: \n" << orb_axes << "\n\n";
    }

    ++i_orbit;
  }
  // std::cout << "axes: \n" << axes << "\n";
  if (axes.cols() == 0) {
    // Strategy 4
    // std::cout << "Strategy 4: (fallback)" << std::endl;
    result = irrep_subspace.colPivHouseholderQr().solve(
        vector_space_prepare(irrep_subspace, vec_compare_tol));
    i_strategy = 4;
  } else {
    result = irrep_subspace.colPivHouseholderQr().solve(axes);
  }

  if (log.has_value()) {
    if (i_strategy == 1) {
      log->indent() << "Success: Strategy 1 (Found orthogonal axes using "
                       "directions from a single orbit)"
                    << std::endl;
      log->indent() << "Orbit=" << i_orb_orthog << ", mult=" << min_mult
                    << std::endl
                    << std::endl;
    } else if (i_strategy == 2) {
      log->indent() << "Success: Strategy 2 (Found orthogonal axes using "
                       "directions from multiple orbits)"
                    << std::endl;
      log->indent() << "Orbit providing each axis: " << i_tot_orthog
                    << std::endl
                    << std::endl;
    }
    // qr decomposition on lowest-multiplicity
    // orbit to find spanning set of axes
    else if (i_strategy == 3) {
      log->indent() << "Success: Strategy 3 (Found orthogonal axes from a QR "
                       "decomposition of the lowest-multiplicity orbit)"
                    << std::endl;
      log->indent() << "Orbit used=" << i_min_mult_orbit
                    << ", mult=" << min_mult << std::endl
                    << std::endl;
    } else if (i_strategy == 4) {
      log->indent() << "Fallback: Strategy 4 (Found orthogonal axes from a QR "
                       "decomposition of irrep subspace matrix)"
                    << std::endl
                    << std::endl;
    }
  }

  // std::cout << "result: \n" << result << "\n"
  //<< "irrep_subspace*result: \n" << irrep_subspace *result << "\n";
  // std::cout << "DONE" << std::endl;
  return result;  //.transpose();
}

}  // namespace irreps

}  // namespace CASM
