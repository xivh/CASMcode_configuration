#ifndef CASM_IntegralCluster_json_io
#define CASM_IntegralCluster_json_io

#include <memory>
#include <optional>

#include "casm/configuration/clusterography/IntegralCluster.hh"

namespace CASM {

namespace clust {
class IntegralCluster;
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

/// \brief Write IntegralCluster to JSON object
jsonParser &to_json(
    clust::IntegralCluster const &clust, jsonParser &json,
    xtal::BasicStructure const &prim,
    std::optional<clust::IntegralCluster> phenomenal = std::nullopt);

/// \brief Read from JSON
void from_json(clust::IntegralCluster &clust, jsonParser const &json,
               xtal::BasicStructure const &prim);

template <>
struct jsonConstructor<clust::IntegralCluster> {
  /// \brief Construct from JSON
  static clust::IntegralCluster from_json(jsonParser const &json,
                                          xtal::BasicStructure const &prim);
};

/// \brief Parse IntegralCluster from JSON
void parse(InputParser<clust::IntegralCluster> &parser,
           xtal::BasicStructure const &prim);

}  // namespace CASM

#endif
