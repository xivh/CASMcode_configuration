#ifndef CASM_irreps_VectorSpaceSymReport
#define CASM_irreps_VectorSpaceSymReport

#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/configuration/irreps/IrrepWedge.hh"

namespace CASM {
namespace irreps {

///\brief Summary of data associated with the action of a symmetry group on a
/// vector space
struct VectorSpaceSymReport {
  /// \brief Constructor
  VectorSpaceSymReport(std::vector<Eigen::MatrixXd> _symgroup_rep,
                       std::vector<IrrepInfo> _irreps,
                       std::vector<SubWedge> _irreducible_wedge,
                       Eigen::MatrixXd const &_symmetry_adapted_subspace,
                       std::vector<std::string> _axis_glossary);

  /// \brief Matrix representation for each operation in the group -- defines
  /// action of group on vector space
  std::vector<Eigen::MatrixXd> symgroup_rep;

  /// \brief A list of all irreducible representation that make up the full
  /// representation
  std::vector<IrrepInfo> irreps;

  /// \brief Names for the irreps
  ///
  /// Form is "irrep_{i}_{j}", where:
  /// - i: irrep group index (starting from 1) distinguishing groups of
  /// identical irreps
  ///   (irreps are identical if they have the same character vectors)
  /// - j: equiv irrep index (starting from 1) distinguishing irreps with the
  /// same
  ///   character vectors
  std::vector<std::string> irrep_names;

  /// \brief Indices of columns in symmetry_adapted_subspace corresponding to
  /// each irrep
  ///
  /// Column c = this->irrep_axes_indices[i][j]: is the column index (starting
  /// with 0) in symmetry_adapted_subspace corresponding to the j-th dim of
  /// this->irreps[i]
  std::vector<std::vector<Index>> irrep_axes_indices;

  /// \brief Irreducible wedge axes in the full vector space by irrep
  ///
  /// irrep_wedge_axes[i] corresponds to the wedge for this->irreps[i]
  std::vector<Eigen::MatrixXd> irrep_wedge_axes;

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
    IrrepDecomposition const &irrep_decomposition, bool calc_wedges = false,
    std::optional<std::vector<std::string>> axis_glossary = std::nullopt);

}  // namespace irreps
}  // namespace CASM

#endif
