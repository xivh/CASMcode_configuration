#ifndef CASM_irreps_IrrepDecomposition_json_io
#define CASM_irreps_IrrepDecomposition_json_io

namespace CASM {

class jsonParser;

namespace irreps {

struct IrrepInfo;

jsonParser &to_json(IrrepInfo const &irrep, jsonParser &json);

}  // namespace irreps
}  // namespace CASM

#endif
