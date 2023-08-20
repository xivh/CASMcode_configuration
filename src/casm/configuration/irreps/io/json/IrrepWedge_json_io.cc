#include "casm/configuration/irreps/io/json/IrrepWedge_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/irreps/IrrepWedge.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace irreps {

jsonParser &to_json(SubWedge const &wedge, jsonParser &json) {
  json["full_wedge_axes"] = wedge.trans_mat.transpose();

  for (Index i = 0; i < wedge.irrep_wedges.size(); ++i) {
    std::string irrep_name =
        "irrep_" + to_sequential_string(i + 1, wedge.irrep_wedges.size());
    json["irrep_wedge_axes"][irrep_name] =
        wedge.irrep_wedges[i].axes.transpose();
  }
  return json;
}

}  // namespace irreps
}  // namespace CASM
