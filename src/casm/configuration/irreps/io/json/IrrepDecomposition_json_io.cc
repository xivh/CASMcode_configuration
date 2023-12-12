#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace irreps {

/// \brief Output irrep characters, if a pseudo irrep.
///
/// If a pseudo irrep, output the "is_direct_sum_of_complex_irreps" portion of
/// the IrrepInfo as JSON.
jsonParser &add_pseudo_irrep_characters(IrrepInfo const &irrep,
                                        jsonParser &json) {
  if (irrep.pseudo_irrep) {
    to_json_array(irrep.characters.imag(),
                  json["is_direct_sum_of_complex_irreps"]
                      ["symop_characters_imag_component"]);
    to_json_array(irrep.characters.real(),
                  json["is_direct_sum_of_complex_irreps"]
                      ["symop_characters_real_component"]);
  }
  return json;
}

/// \brief Represent IrrepInfo as JSON
///
/// \param irrep IrrepInfo to represent as JSON
/// \param json Expects an existing JSON object
/// \return Reference to modified JSON object
///
/// \code
/// {
///   "is_direct_sum_of_complex_irreps": optional, object
///     Only included for pseudo irreps. If irrep is real but was created as
///     direct sum of two complex irreps in this case, the 'irrep' is reducible,
///     this pseudo irrep is the most-reduced representation that can still
///     have real basis vectors.
///
///     "symop_characters_imag_component": optional, array of float
///       The imaginary component of the irrep characters.
///
///      "symop_characters_real_component": optional, array of float
///       The real component of the irrep characters.
///
///   "axes": object
///       An `irrep_dim` x `vector_dim` matrix transforms a vector from the
///       initial vector space into a vector in the irreducible vector space.
///       The transformation matrix may be complex, so the real and imaginary
///       components are included separately.
///
///     "real": array_like of float
///       Real component of the irrep transformation matrix.
///
///     "imag": optional, array_like of float
///       Negative of the imaginary component of the irrep transformation
///       matrix, if it is complex.
///
///   "high_symmetry_directions": optional, array_like
///     Vectors in the irreducible vector space that correspond to high-symmetry
///     directions. X[i] is the i'th orbit of equivalent high-symmetry
///     directions and X[i].size() is the symmetric multiplicity of a direction
///     in that orbit. X[i][j] is an array of float, specifying a high-symmetry
///     direction in the irreducible vector space.
///
//
/// \endcode
jsonParser &to_json(IrrepInfo const &irrep, jsonParser &json) {
  add_pseudo_irrep_characters(irrep, json);

  if (!almost_zero(irrep.trans_mat.imag())) {
    json["axes"]["imaginary"] = -irrep.trans_mat.imag();
  }

  json["axes"]["real"] = irrep.trans_mat.real();

  if (!irrep.directions.empty()) {
    json["high_symmetry_directions"].put_array(irrep.directions.size());
  }

  for (Index i = 0; i < irrep.directions.size(); ++i) {
    json["high_symmetry_directions"][i].put_array(irrep.directions[i].size());
    for (Index j = 0; j < irrep.directions[i].size(); ++j) {
      to_json_array(
          Eigen::MatrixXd(irrep.trans_mat.real() * irrep.directions[i][j]),
          json["high_symmetry_directions"][i][j]);
    }
  }
  return json;
}

}  // namespace irreps
}  // namespace CASM
