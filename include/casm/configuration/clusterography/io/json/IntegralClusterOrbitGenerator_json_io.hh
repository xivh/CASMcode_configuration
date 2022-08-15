#ifndef CASM_clust_IntegralClusterOrbitGenerator_json_io
#define CASM_clust_IntegralClusterOrbitGenerator_json_io

#include <vector>

namespace CASM {

namespace clust {
struct IntegralClusterOrbitGenerator;
}  // namespace clust

namespace xtal {
class BasicStructure;
}

template <typename T>
class InputParser;
template <typename T>
struct jsonConstructor;
class jsonParser;

/// Write custom orbit specs to JSON
jsonParser &to_json(const clust::IntegralClusterOrbitGenerator &orbit_generator,
                    jsonParser &json, xtal::BasicStructure const &prim);

/// Parse vector of IntegralClusterOrbitGenerator ("orbit_specs") from JSON
void parse(
    InputParser<std::vector<clust::IntegralClusterOrbitGenerator>> &parser,
    xtal::BasicStructure const &prim);

}  // namespace CASM

#endif
