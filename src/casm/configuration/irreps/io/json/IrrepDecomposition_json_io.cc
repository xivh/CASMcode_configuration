#include "casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/IrrepDecomposition.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace irreps {

jsonParser &to_json(IrrepInfo const &irrep, jsonParser &json) {
  if (irrep.pseudo_irrep) {
    to_json_array(irrep.characters.imag(),
                  json["is_direct_sum_of_complex_irreps"]
                      ["symop_characters_imag_component"]);
    to_json_array(irrep.characters.real(),
                  json["is_direct_sum_of_complex_irreps"]
                      ["symop_characters_real_component"]);
  }

  /*
  if(!almost_zero(irrep.trans_mat.imag())) {
    json["axes"]["imaginary"] = -irrep.trans_mat.imag();
  }

  json["axes"]["real"] = irrep.trans_mat.real();
  */
  /*
  if(irrep.directions.empty()) {
    json["high_symmetry_directions"] = "none";
  }
  else {
    json["high_symmetry_directions"].put_array(irrep.directions.size());
  }

  for(Index i = 0; i < irrep.directions.size(); ++i) {
    json["high_symmetry_directions"][i].put_array(irrep.directions[i].size());
    for(Index j = 0; j < irrep.directions[i].size(); ++j) {
      to_json_array(Eigen::MatrixXd(irrep.trans_mat.real()*irrep.directions[i][j]),
  json["high_symmetry_directions"][i][j]);
    }
  }
  */
  return json;
}

}  // namespace irreps
}  // namespace CASM
