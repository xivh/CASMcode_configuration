#include "casm/configuration/irreps/VectorSpaceSymReport.hh"

namespace CASM {
namespace irreps {

/// Construct VectorSpaceSymReport
///
/// \param irrep_decomposition An IrrepDecomposition
/// \param calc_wedges If true, 'irreducible_wedge' is constructed. If false,
///     'irreducible_wedge' is empty.
VectorSpaceSymReport vector_space_sym_report(
    IrrepDecomposition const &irrep_decomposition, bool calc_wedges) {
  VectorSpaceSymReport result;
  result.irreps = irrep_decomposition.irreps;
  result.symmetry_adapted_subspace =
      irrep_decomposition.symmetry_adapted_subspace;
  for (Index element_index : irrep_decomposition.head_group) {
    auto matrix_rep = irrep_decomposition.fullspace_rep[element_index];
    result.symgroup_rep.push_back(matrix_rep);
  }
  result.axis_glossary =
      std::vector<std::string>(result.symmetry_adapted_subspace.rows(), "x");
  Index i = 0;
  for (std::string &x : result.axis_glossary) {
    x += std::to_string(++i);
  }
  if (calc_wedges) {
    result.irreducible_wedge = make_symrep_subwedges(irrep_decomposition);
  }
  return result;
}

}  // namespace irreps
}  // namespace CASM
