#ifndef CASM_occ_events_OccEventCounter_json_io
#define CASM_occ_events_OccEventCounter_json_io

namespace CASM {
namespace occ_events {
struct OccEventCounterStateInfo;
struct OccSystem;
}  // namespace occ_events

class jsonParser;

jsonParser &to_json(occ_events::OccEventCounterStateInfo const &state_info,
                    jsonParser &json, occ_events::OccSystem const &system);

}  // namespace CASM

#endif
