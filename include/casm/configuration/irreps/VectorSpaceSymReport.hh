#ifndef CASM_irreps_VectorSpaceSymReport
#define CASM_irreps_VectorSpaceSymReport

#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/configuration/irreps/IrrepWedge.hh"

namespace CASM {
namespace irreps {

///\brief Summary of data associated with the action of a symmetry group on a
/// vector space
struct VectorSpaceSymReport {
  /// \brief Matrix representation for each operation in the group -- defines
  /// action of group on vector space
  std::vector<Eigen::MatrixXd> symgroup_rep;

  /// \brief A list of all irreducible representation that make up the full
  /// representation
  std::vector<IrrepInfo> irreps;

  /// \brief Irreducible wedge in the vector space
  /// encoded as a vector of symmetrically distinct SubWedges
  std::vector<SubWedge> irreducible_wedge;

  /// \brief Symmetry-oriented subspace of the vector space (columns are the
  /// basis vectors)
  Eigen::MatrixXd symmetry_adapted_subspace;

  /// \brief Names given to individual axes in initial (un-adapted) vector
  /// space, corresponding to rows of symmetry_adapted_dof_subspace
  std::vector<std::string> axis_glossary;
};

/// Construct VectorSpaceSymReport
VectorSpaceSymReport vector_space_sym_report(
    IrrepDecomposition const &irrep_decomposition, bool calc_wedges = false);

}  // namespace irreps
}  // namespace CASM

#endif
