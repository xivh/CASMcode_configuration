#include "casm/configuration/clusterography/io/json/EquivalentsInfo_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"

namespace CASM {

/// Complete "equivalents_info.json" format:
///
/// "equivalent_generating_ops": array of int (required)
///     Indices of factor group operations used to generate
///     equivalent phenomenal clusters
/// "equivalents": array of LocalOrbitTree
///     An array of local orbit trees
/// "factor_group": SymGroup
///     The prim factor group
/// "prototype": LocalOrbitTree
///     The prototype local orbit tree
///
///
/// LocalOrbitTree format:
///
/// "orbits": array of Cluster Orbit
///     Array of local cluster orbits
/// "phenomenal": IntegralCluster
///     The phenomenal cluster of the local orbits
/// "prim": Prim
///     The primitive crystal structure and degrees of freedom (DoF)
///
///
/// Cluster Orbit format:
///
/// "linear_orbit_index": int
///     Linear index of orbit in array of orbits
/// "mult": int
///     Number of elements in the orbit
/// "prototype": IntegralCluster
///     The prototype cluster
/// "elements": list[IntegralCluster]
///     Array of cluster
///

/// \brief Output minimal "equivalents info" to JSON
jsonParser &to_json(clust::EquivalentsInfo const &equivalents_info,
                    jsonParser &json, xtal::BasicStructure const &prim) {
  json["equivalent_generating_ops"] =
      equivalents_info.equivalent_generating_op_indices;
  json["equivalents"] = jsonParser::array();
  for (auto const &phenomenal_cluster : equivalents_info.phenomenal_clusters) {
    jsonParser equivalent;
    to_json(phenomenal_cluster, equivalent["phenomenal"], prim);
    json["equivalents"].push_back(equivalent);
  }
  return json;
}

/// \brief Parse minimal clexulator::EquivalentsInfo from JSON
void parse(InputParser<clust::EquivalentsInfo> &parser,
           xtal::BasicStructure const &prim) {
  parser.value = std::make_unique<clust::EquivalentsInfo>();
  auto &equivalents_info = *parser.value;

  parser.require(equivalents_info.equivalent_generating_op_indices,
                 "equivalent_generating_ops");

  if (parser.self.find("equivalents") == parser.self.end() ||
      !parser.self["equivalents"].is_array()) {
    parser.error.insert("Error: missing 'equivalents' array ");
    parser.value = nullptr;
    return;
  }

  clust::IntegralCluster cluster;

  auto const &equivalent_json = parser.self["equivalents"];
  for (Index i = 0; i < equivalent_json.size(); ++i) {
    parser.require(cluster,
                   fs::path("equivalents") / std::to_string(i) / "phenomenal",
                   prim);
    equivalents_info.phenomenal_clusters.push_back(cluster);
  }

  if (!parser.valid()) {
    parser.value = nullptr;
  }
}

}  // namespace CASM
