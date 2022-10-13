#ifndef CASM_occ_events_OccSystem_json_io
#define CASM_occ_events_OccSystem_json_io

#include <memory>

namespace CASM {
namespace occ_events {
struct OccSystem;
}
namespace xtal {
class BasicStructure;
}

template <typename T>
class InputParser;
class jsonParser;

jsonParser &to_json(occ_events::OccSystem const &system, jsonParser &json);

void parse(InputParser<occ_events::OccSystem> &parser,
           std::shared_ptr<xtal::BasicStructure const> const &prim);

}  // namespace CASM

#endif
