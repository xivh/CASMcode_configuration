#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {

/// \brief Output irrep characters, if a pseudo irrep.
///
/// If a pseudo irrep, output the "is_direct_sum_of_complex_irreps" portion of
/// the IrrepInfo as JSON.
jsonParser &add_pseudo_irrep_characters(irreps::IrrepInfo const &irrep,
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
///   "irrep_type": optional, integer
///     Index that is the same for irreps with approximately the same
///     characters.
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
jsonParser &to_json(irreps::IrrepInfo const &irrep, jsonParser &json) {
  // Transformation matrix / axes
  if (!almost_zero(irrep.trans_mat.imag())) {
    json["axes"]["imaginary"] = -irrep.trans_mat.imag();
  }
  json["axes"]["real"] = irrep.trans_mat.real();

  // Characters
  if (!almost_zero(irrep.characters.imag())) {
    json["characters"]["imaginary"] = irrep.characters.imag();
  }
  json["characters"]["real"] = irrep.characters.real();

  json["complex"] = irrep.complex;

  // Pseudo irrep characters
  add_pseudo_irrep_characters(irrep, json);
  json["pseudo_irrep"] = irrep.pseudo_irrep;

  // Index
  if (irrep.index.has_value()) {
    json["index"] = *irrep.index;
  }

  // Irrep type
  if (irrep.irrep_type.has_value()) {
    json["irrep_type"] = *irrep.irrep_type;
  }

  // Frobenius-Schur indicator
  json["frobenius_schur_indicator"] = irrep.frobenius_schur_indicator;

  // High-symmetry directions
  if (irrep.directions.has_value() && !irrep.directions->empty()) {
    json["high_symmetry_directions"].put_array(irrep.directions->size());
    for (Index i = 0; i < irrep.directions->size(); ++i) {
      json["high_symmetry_directions"][i].put_array(
          (*irrep.directions)[i].size());
      for (Index j = 0; j < (*irrep.directions)[i].size(); ++j) {
        to_json_array(
            Eigen::MatrixXd(irrep.trans_mat.real() * (*irrep.directions)[i][j]),
            json["high_symmetry_directions"][i][j]);
      }
    }
  }
  return json;
}

/// Read irreps::IrrepInfo from JSON
irreps::IrrepInfo jsonConstructor<irreps::IrrepInfo>::from_json(
    jsonParser const &json) {
  InputParser<irreps::IrrepInfo> parser(json);
  std::stringstream ss;
  ss << "Error: Invalid IrrepInfo JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Parse IrrepInfo from JSON with error messages
void parse(InputParser<irreps::IrrepInfo> &parser) {
  // Read axes / transformation matrix
  Eigen::MatrixXd trans_mat_real;
  fs::path axes_path{"axes"};
  parser.require(trans_mat_real, axes_path / "real");

  Eigen::MatrixXd trans_mat_imag =
      Eigen::MatrixXd::Zero(trans_mat_real.rows(), trans_mat_real.cols());
  fs::path imag_path = axes_path / "imaginary";
  parser.optional(trans_mat_imag, imag_path);

  Eigen::MatrixXcd trans_mat(trans_mat_real.rows(), trans_mat_real.cols());
  trans_mat.real() = trans_mat_real;
  trans_mat.imag() = trans_mat_imag;

  // Read characters / pseudo irrep info
  bool pseudo_irrep;
  parser.require(pseudo_irrep, "pseudo_irrep");

  fs::path char_path{"characters"};
  Eigen::VectorXd char_real;
  parser.require(char_real, char_path / "real");

  Eigen::VectorXd char_imag = Eigen::VectorXd::Zero(char_real.size());
  fs::path char_imag_path = char_path / "imaginary";
  parser.optional(char_imag, char_imag_path);

  Eigen::VectorXcd characters(char_real.size());
  characters.real() = char_real;
  characters.imag() = char_imag;

  // Index
  std::optional<Index> index = std::nullopt;
  if (parser.self.find("index") != parser.self.end()) {
    Index index_val;
    parser.require(index_val, fs::path{"index"});
    index = index_val;
  }

  // Irrep type
  std::optional<Index> irrep_type = std::nullopt;
  if (parser.self.find("irrep_type") != parser.self.end()) {
    Index irrep_type_val;
    parser.require(irrep_type_val, fs::path{"irrep_type"});
    irrep_type = irrep_type_val;
  }

  // High-symmetry directions
  std::optional<std::vector<std::vector<Eigen::VectorXd>>> directions =
      std::nullopt;
  fs::path directions_path{"high_symmetry_directions"};
  if (parser.self.find(directions_path) != parser.self.end()) {
    jsonParser directions_json = parser.self[directions_path];
    directions = std::vector<std::vector<Eigen::VectorXd>>{};
    directions->resize(directions_json.size());
    for (Index i = 0; i < directions_json.size(); ++i) {
      jsonParser orbit_json = directions_json[i];
      (*directions)[i].resize(orbit_json.size());
      for (Index j = 0; j < orbit_json.size(); ++j) {
        Eigen::VectorXd direction;
        parser.require(direction,
                       directions_path / std::to_string(i) / std::to_string(j));
        (*directions)[i][j] = direction;
      }
    }
  }

  // Frobenius-Schur indicator
  int frobenius_schur_indicator;
  parser.require(frobenius_schur_indicator, "frobenius_schur_indicator");

  if (parser.valid()) {
    parser.value =
        notstd::make_unique<irreps::IrrepInfo>(trans_mat, characters);
    parser.value->pseudo_irrep = pseudo_irrep;
    parser.value->index = index;
    parser.value->irrep_type = irrep_type;
    parser.value->directions = directions;
    parser.value->frobenius_schur_indicator = frobenius_schur_indicator;
  }
}

}  // namespace CASM
