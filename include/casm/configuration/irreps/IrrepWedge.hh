#ifndef CASM_irreps_IrrepWedge
#define CASM_irreps_IrrepWedge

#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/configuration/irreps/definitions.hh"

namespace CASM {
namespace irreps {

/// An irreducible wedge in an irreducible vector space
struct IrrepWedge {
  IrrepWedge(IrrepInfo _irrep_info, Eigen::MatrixXd _axes);

  /// The description of the associated irreducible vector space
  IrrepInfo irrep_info;

  /// columns of 'axes' are high-symmetry direction that form the edges
  /// of the irreducible wedge. The number of columns equals the
  /// dimension of the irreducible subspace. The number of rows is usually the
  /// fullspace dimension, but depending on context it could be a subspace.
  Eigen::MatrixXd axes;

  /// Symmetric multiplicity of each column of 'axes'
  std::vector<Index> mult;
};

/// Construct a "dummy" IrrepWedge with user specified axes
IrrepWedge make_dummy_irrep_wedge(Eigen::MatrixXd const &axes);

/// SubWedge is a vector of IrrepWedge, one from each irreducible subspace
///
/// Together, the IrrepWedges that comprise the Subwedge span the entire space
/// However, it is likely that the orbit of equivalent SubWedges does not *fill*
/// the entire space. Then multiple SubWedge combine to fill the entire space.
class SubWedge {
 public:
  SubWedge(std::vector<IrrepWedge> const &_irrep_wedges);

  /// IrrepWedges comprising the Subwedge
  std::vector<IrrepWedge> irrep_wedges;

  /// Transformation matrix to convert from a vector in terms of the SubWedge
  /// axes to a vector in the original vector space
  Eigen::MatrixXd trans_mat;

 private:
  static Eigen::MatrixXd _subwedge_to_trans_mat(
      std::vector<IrrepWedge> const &_irrep_wedges);
};

/// Makes a "dummy" SubWedge from a single "dummy" IrrepWedge with given axes
SubWedge make_dummy_subwedge(Eigen::MatrixXd const &axes);

/// Make IrrepWedges from an IrrepDecomposition
std::vector<IrrepWedge> make_irrep_wedges(
    IrrepDecomposition const &irrep_decomposition);

/// \brief Find full irreducible wedge of a group-represented vector space, as
/// a vector of SubWedges, from an IrrepDecomposition
std::vector<SubWedge> make_symrep_subwedges(
    IrrepDecomposition const &irrep_decomposition);

}  // namespace irreps
}  // namespace CASM

#endif
