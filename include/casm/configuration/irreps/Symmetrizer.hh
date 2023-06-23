#ifndef CASM_irreps_Symmetrizer
#define CASM_irreps_Symmetrizer

#include "casm/configuration/irreps/definitions.hh"

namespace CASM {
namespace irreps {

/// Find high-symmetry directions in a irreducible space
multivector<Eigen::VectorXcd>::X<2> make_irrep_special_directions(
    MatrixRep const &rep, GroupIndices const &head_group,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol,
    GroupIndicesOrbitVector const &cyclic_subgroups,
    GroupIndicesOrbitVector const &all_subgroups,
    bool use_all_subgroups = false);

/// Make an irreducible space symmetrizer matrix using special directions
Eigen::MatrixXcd make_irrep_symmetrizer_matrix(
    multivector<Eigen::VectorXcd>::X<2> const &irrep_special_directions,
    Eigen::MatrixXcd const &irrep_subspace, double vec_compare_tol);

}  // namespace irreps
}  // namespace CASM

#endif
