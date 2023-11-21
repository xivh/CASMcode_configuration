#include "casm/configuration/irreps/VectorSpaceSymReport.hh"

#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace irreps {

/// \brief Constructor
VectorSpaceSymReport::VectorSpaceSymReport(
    std::vector<Eigen::MatrixXd> _symgroup_rep, std::vector<IrrepInfo> _irreps,
    std::vector<SubWedge> _irreducible_wedge,
    Eigen::MatrixXd const &_symmetry_adapted_subspace,
    std::vector<std::string> _axis_glossary)
    : symgroup_rep(std::move(_symgroup_rep)),
      irreps(std::move(_irreps)),
      irreducible_wedge(std::move(_irreducible_wedge)),
      symmetry_adapted_subspace(_symmetry_adapted_subspace),
      axis_glossary(std::move(_axis_glossary)) {
  // ~~~ Identify irrep_names, irrep_axes_indices, irrep_wedges  ~~~
  std::vector<Index> mults;
  for (auto const &irrep : this->irreps) {
    if (irrep.index == 0) mults.push_back(0);
    mults.back()++;
  }

  Index NQ = this->symmetry_adapted_subspace.cols();
  Index i = 0;  // irreps with unique characters index
  Index l = 0;  // irrep index
  Index q = 0;  // symmetry_adapted_subspace column index
  for (Index mult : mults) {
    ++i;
    for (Index m = 0; m < mult; ++m) {  // irreps with same character index

      // irrep_names
      std::string irrep_name = "irrep_" +
                               to_sequential_string(i, mults.size()) + "_" +
                               to_sequential_string(m + 1, mult);
      this->irrep_names.push_back(irrep_name);

      // irrep_wedge_axes
      auto const &irrep = this->irreps[l];
      if (this->irreducible_wedge.size()) {
        this->irrep_wedge_axes.push_back(
            this->irreducible_wedge[0].irrep_wedges[l].axes);
      }

      // irrep_axes_indices
      std::vector<Index> irrep_axes_indices;
      for (Index a = 0; a < irrep.irrep_dim; ++a, ++q) {
        irrep_axes_indices.push_back(q);
      }
      this->irrep_axes_indices.push_back(irrep_axes_indices);

      // increment irrep index
      ++l;
    }
  }
}

/// Construct VectorSpaceSymReport
///
/// \param irrep_decomposition An IrrepDecomposition
/// \param calc_wedges If true, 'irreducible_wedge' is constructed. If false,
///     'irreducible_wedge' is empty.
/// \param axis_glossary If has value, copied to
/// VectorSpaceSymReport.axis_glossary;
///     otherwise, axis_glossary is set to {"x1", "x2", ...}
VectorSpaceSymReport vector_space_sym_report(
    IrrepDecomposition const &irrep_decomposition, bool calc_wedges,
    std::optional<std::vector<std::string>> axis_glossary) {
  std::vector<Eigen::MatrixXd> symgroup_rep;
  for (Index element_index : irrep_decomposition.head_group) {
    symgroup_rep.push_back(irrep_decomposition.fullspace_rep[element_index]);
  }

  if (!axis_glossary.has_value()) {
    axis_glossary = std::vector<std::string>(
        irrep_decomposition.symmetry_adapted_subspace.rows(), "x");
    Index i = 0;
    for (std::string &x : axis_glossary.value()) {
      x += std::to_string(++i);
    }
  }

  std::vector<SubWedge> irreducible_wedge;
  if (calc_wedges) {
    irreducible_wedge = make_symrep_subwedges(irrep_decomposition);
  }

  return VectorSpaceSymReport(
      symgroup_rep, irrep_decomposition.irreps, irreducible_wedge,
      irrep_decomposition.symmetry_adapted_subspace, axis_glossary.value());
}

}  // namespace irreps
}  // namespace CASM
