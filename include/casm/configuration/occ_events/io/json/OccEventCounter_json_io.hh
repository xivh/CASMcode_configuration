#ifndef CASM_occ_events_OccEventCounter_json_io
#define CASM_occ_events_OccEventCounter_json_io

namespace CASM {
namespace occ_events {
struct OccEventCounterParameters;
struct OccEventCounterStateInfo;
struct OccSystem;
}  // namespace occ_events

template <typename T>
class InputParser;
class jsonParser;

jsonParser &to_json(occ_events::OccEventCounterStateInfo const &state_info,
                    jsonParser &json, occ_events::OccSystem const &system);

jsonParser &to_json(occ_events::OccEventCounterParameters const &params,
                    jsonParser &json);

void parse(InputParser<occ_events::OccEventCounterParameters> &parser);

}  // namespace CASM

#endif
