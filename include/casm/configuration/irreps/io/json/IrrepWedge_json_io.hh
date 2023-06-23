#ifndef CASM_irreps_IrrepWedge_json_io
#define CASM_irreps_IrrepWedge_json_io

namespace CASM {

class jsonParser;

namespace irreps {

class SubWedge;

jsonParser &to_json(SubWedge const &wedge, jsonParser &json);

}  // namespace irreps
}  // namespace CASM

#endif
