#include "casm/configuration/occ_events/io/json/OccEventCounter_json_io.hh"

#include <optional>

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"

namespace CASM {

jsonParser &to_json(occ_events::OccEventCounterStateInfo const &state_info,
                    jsonParser &json, occ_events::OccSystem const &system) {
  to_json(state_info.occ_init, json["occ_init"]);
  to_json(state_info.occ_final, json["occ_final"]);
  to_json(state_info.position_init, json["position_init"], system);
  to_json(state_info.position_final, json["position_final"], system);
  to_json(state_info.cluster, json["cluster"], *system.prim);
  to_json(state_info.occ_event, json["event"], system);
  if (state_info.fails.size()) {
    to_json(state_info.fails, json["fails"]);
  } else {
    to_json(true, json["passes"]);
  }
  return json;
}

}  // namespace CASM
