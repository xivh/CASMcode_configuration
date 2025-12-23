#ifndef CASM_irreps_VectorSpaceSymReport_json_io
#define CASM_irreps_VectorSpaceSymReport_json_io

namespace CASM {

class jsonParser;

namespace irreps {

struct VectorSpaceSymReport;

jsonParser &to_json(VectorSpaceSymReport const &obj, jsonParser &json,
                    bool include_symop_matrices = true);

}  // namespace irreps
}  // namespace CASM

#endif
