#ifndef CASM_clust_EquivalentsInfo_json_io
#define CASM_clust_EquivalentsInfo_json_io

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/orbits.hh"

namespace CASM {

template <typename T>
class InputParser;
class jsonParser;

/// \brief Output minimal "equivalents info" to JSON
jsonParser &to_json(clust::EquivalentsInfo const &equivalents_info,
                    jsonParser &json, xtal::BasicStructure const &prim);

/// \brief Parse minimal clexulator::EquivalentsInfo from JSON
void parse(InputParser<clust::EquivalentsInfo> &parser,
           xtal::BasicStructure const &prim);

}  // namespace CASM

#endif
