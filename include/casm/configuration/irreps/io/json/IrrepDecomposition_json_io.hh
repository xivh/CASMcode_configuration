#ifndef CASM_irreps_IrrepDecomposition_json_io
#define CASM_irreps_IrrepDecomposition_json_io

namespace CASM {

template <typename T>
struct jsonConstructor;
class jsonParser;
template <typename T>
class InputParser;

namespace irreps {

struct IrrepInfo;
struct IrrepDecomposition;

}  // namespace irreps

/// \brief Output irrep characters, if a pseudo irrep.
jsonParser &add_pseudo_irrep_characters(irreps::IrrepInfo const &irrep,
                                        jsonParser &json);

/// \brief Represent IrrepInfo as JSON
jsonParser &to_json(irreps::IrrepInfo const &irrep, jsonParser &json);

template <>
struct jsonConstructor<irreps::IrrepInfo> {
  /// Read irreps::IrrepInfo from JSON
  static irreps::IrrepInfo from_json(jsonParser const &json);
};

/// \brief Parse IrrepInfo from JSON with error messages
void parse(InputParser<irreps::IrrepInfo> &parser);

}  // namespace CASM

#endif
