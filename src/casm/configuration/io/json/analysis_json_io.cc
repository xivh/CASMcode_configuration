#include "casm/configuration/io/json/analysis_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clexulator/io/json/DoFSpace_json_io.hh"
#include "casm/configuration/config_space_analysis.hh"
#include "casm/configuration/dof_space_analysis.hh"
#include "casm/configuration/io/json/Configuration_json_io.hh"
#include "casm/configuration/irreps/io/json/VectorSpaceSymReport_json_io.hh"

namespace CASM {
/// \brief Write ConfigSpaceAnalysisResults to JSON
jsonParser &to_json(config::ConfigSpaceAnalysisResults const &results,
                    jsonParser &json) {
  json.put_obj();
  to_json(results.standard_dof_space, json["standard_dof_space"]);
  to_json(results.equivalent_dof_values, json["equivalent_dof_values"],
          jsonParser::as_array());
  to_json(results.equivalent_configurations, json["equivalent_configurations"]);
  to_json(results.projector, json["projector"]);
  to_json(results.eigenvalues, json["eigenvalues"], jsonParser::as_array());
  to_json(results.symmetry_adapted_dof_space,
          json["symmetry_adapted_dof_space"]);
  return json;
}

/// \brief Write DoFSpaceAnalysisResults to JSON
jsonParser &to_json(config::DoFSpaceAnalysisResults const &results,
                    jsonParser &json) {
  json.put_obj();
  to_json(results.symmetry_report, json["symmetry_report"]);
  to_json(results.symmetry_adapted_dof_space,
          json["symmetry_adapted_dof_space"]);
  return json;
}

}  // namespace CASM
