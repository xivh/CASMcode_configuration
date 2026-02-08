#include "casm/configuration/irreps/Symmetrizer.hh"

#include <vector>

#include "casm/casm_io/container/stream_io.hh"
#include "casm/configuration/irreps/SimpleOrbit_impl.hh"
#include "casm/configuration/irreps/VectorSymCompare_v2.hh"
#include "casm/global/threads.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

namespace irreps {

/// Find high-symmetry directions in a irreducible space
///
/// Notes:
/// - This method is multithreaded, parallelizing the loop over subgroup orbits
///
/// \param rep Matrix representation of head_group, this defines group action
/// on the underlying vector space
/// \param head_group Group for which the irreps are to be found
/// \param irrep_subspace A column vector matrix representing a basis of the
///     irreducible space in which high symmetry directions will be found. The
///     number of rows must equal `rep.dim()`, the number of columns is equal to
///     the dimension of the irreducible space.
/// \param vec_compare_tol Tolerance for elementwise floating-point comparisons
///     of vectors
/// \param subgroup_orbits Orbits of subgroups of `head_group`, where each
///     subgroup is represented as a set of indices into `rep`. The vector
///     `subgroup_orbits[i][j]` contains the indices of elements in the `j`-th
///     subgroup of the `i`-th orbit of equivalent subgroups.
/// \param log Optional Log object for logging progress
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
    GroupIndicesOrbitSet const &subgroup_orbits, std::optional<Log> log) {
  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Get subgroup indices...";
    append_time(*log, 1);
  }

  GroupIndicesOrbitSet const &sgroups = subgroup_orbits;

  // Copy sgroups to a vector for indexed access in parallel tasks
  std::vector<GroupIndicesOrbit> sgroups_vec;
  for (const auto &orbit : sgroups) {
    sgroups_vec.push_back(orbit);
  }

  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Get subgroup indices: DONE" << std::endl << std::endl;
    log->indent() << "Number of subgroup orbits = " << sgroups.size();
    append_time(*log, 2);
    log->indent() << "Applying Reynolds operator to find special directions... "
                  << std::endl
                  << std::endl;
    log->indent() << "Special directions: ";
  }

  Index dim = rep[0].rows();

  std::set<SimpleOrbit<VectorSymCompare>> orbit_result;
  std::vector<std::set<SimpleOrbit<VectorSymCompare>>> per_thread_orbit_result;
  per_thread_orbit_result.resize(max_threads());

  // Define the worker that finds special directions for a chunk of subgroups
  auto worker = [&](Index start, Index end, Index thread_id) {
    std::set<SimpleOrbit<VectorSymCompare>> &local_orbit_result =
        per_thread_orbit_result[thread_id];
    VectorSymCompare local_sym_compare{rep, vec_compare_tol};
    Eigen::MatrixXcd local_irrep_subspace = irrep_subspace;

    // Local Reynolds matrix
    Eigen::MatrixXd R_local(dim, dim);

    for (Index orbit_idx = start; orbit_idx < end; ++orbit_idx) {
      const auto &orbit = sgroups_vec[orbit_idx];

      // Build Reynolds operator for this subgroup
      R_local.setZero(dim, dim);
      for (Index element_index : *(orbit.begin())) {
        R_local += rep[element_index];
      }

      // Apply Reynolds operator to irrep subspace
      Eigen::MatrixXcd projected = R_local * local_irrep_subspace;

      // If projection is (near) zero, skip
      if (projected.norm() < TOL) return;

      // Find spanning vectors of column space of
      // R*irrep_space
      auto QR = projected.colPivHouseholderQr();
      QR.setThreshold(TOL);

      // If more than one spanning vector, not a
      // unique special direction
      if (QR.rank() > 1) continue;

      Eigen::MatrixXcd Q = QR.matrixQ();
      Eigen::VectorXcd v = Q.col(0);

      local_orbit_result.emplace(v, head_group.begin(), head_group.end(),
                                 local_sym_compare);

      local_orbit_result.emplace(-v, head_group.begin(), head_group.end(),
                                 local_sym_compare);
    }
  };

  threaded_run(sgroups_vec.size(), worker);

  // Insert local_orbit_result elements into shared orbit_result
  for (const auto &orbits : per_thread_orbit_result) {
    for (const auto &orbit : orbits) {
      orbit_result.insert(orbit);
      auto res1 = orbit_result.insert(orbit);
      if (res1.second == true) {
        if (log.has_value() && log->verbosity() >= Log::verbose &&
            log->print()) {
          log->ostream() << "*";
        }
      }

      auto res2 = orbit_result.insert(orbit);
      if (res2.second == true) {
        if (log.has_value() && log->verbosity() >= Log::verbose &&
            log->print()) {
          log->ostream() << "*";
        }
      }
    }
  }

  if (log.has_value() && log->verbosity() >= Log::verbose && log->print()) {
    log->ostream() << std::endl << std::endl;
  }

  multivector<Eigen::VectorXcd>::X<2> result_mt;
  for (auto const &orbit : orbit_result) {
    result_mt.emplace_back(orbit.begin(), orbit.end());
  }

  if (log.has_value() && log->verbosity() >= Log::verbose) {
    log->indent() << "Found " << result_mt.size()
                  << " orbits of special directions." << std::endl
                  << std::endl;
  }

  return result_mt;
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
  if (irrep_special_directions.size() && log.has_value() &&
      log->verbosity() >= Log::verbose) {
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

  if (log.has_value() && log->verbosity() >= Log::verbose) {
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
