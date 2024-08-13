#ifndef CASM_config_analysis_json_io
#define CASM_config_analysis_json_io

namespace CASM {

class jsonParser;

namespace config {
struct ConfigSpaceAnalysisResults;
struct DoFSpaceAnalysisResults;
}  // namespace config

/// \brief Write ConfigSpaceAnalysisResults to JSON
jsonParser &to_json(config::ConfigSpaceAnalysisResults const &results,
                    jsonParser &json);

/// \brief Write DoFSpaceAnalysisResults to JSON
jsonParser &to_json(config::DoFSpaceAnalysisResults const &results,
                    jsonParser &json);

}  // namespace CASM

#endif
