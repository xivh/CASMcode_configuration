#ifndef CASM_occ_events_OccSystem_json_io
#define CASM_occ_events_OccSystem_json_io

namespace CASM {
namespace occ_events {
struct OccSystem;
}

class jsonParser;

jsonParser &to_json(occ_events::OccSystem const &invariants, jsonParser &json);

}  // namespace CASM

#endif
