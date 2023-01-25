#ifndef CASM_clust_ClusterSpecs_json_io
#define CASM_clust_ClusterSpecs_json_io

#include <memory>

#include "casm/configuration/clusterography/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;
template <typename T>
struct jsonConstructor;
class jsonParser;

/// \brief Write ClusterSpecs to JSON object
jsonParser &to_json(clust::ClusterSpecs const &cluster_specs, jsonParser &json);

/// \brief Read from JSON
void from_json(
    clust::ClusterSpecs &cluster_specs, jsonParser const &json,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

template <>
struct jsonConstructor<clust::ClusterSpecs> {
  /// \brief Construct from JSON
  static clust::ClusterSpecs from_json(
      jsonParser const &json,
      std::shared_ptr<xtal::BasicStructure const> const &prim,
      std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
      std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);
};

/// \brief Parse ClusterSpecs from JSON
void parse(
    InputParser<clust::ClusterSpecs> &parser,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep);

}  // namespace CASM

#endif
