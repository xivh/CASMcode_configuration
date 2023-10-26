#ifndef CASM_irreps_IrrepDecomposition_json_io
#define CASM_irreps_IrrepDecomposition_json_io

namespace CASM {

class jsonParser;

namespace irreps {

struct IrrepInfo;

/// \brief Output irrep characters, if a pseudo irrep.
jsonParser &add_pseudo_irrep_characters(IrrepInfo const &irrep,
                                        jsonParser &json);

/// \brief Represent IrrepInfo as JSON
jsonParser &to_json(IrrepInfo const &irrep, jsonParser &json);

}  // namespace irreps
}  // namespace CASM

#endif
