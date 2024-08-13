#ifndef CASM_occ_events_OccEvent_json_io
#define CASM_occ_events_OccEvent_json_io

#include <functional>
#include <memory>
#include <optional>

#include "casm/configuration/occ_events/definitions.hh"
#include "casm/crystallography/io/SymInfo_json_io.hh"
#include "casm/crystallography/io/SymInfo_stream_io.hh"

namespace CASM {

class jsonParser;

template <typename T>
class InputParser;
template <typename T>
struct jsonConstructor;

namespace occ_events {

struct OccEventOutputOptions {
  // --- OccEvent printing options ---
  bool include_cluster = true;
  bool include_cluster_occupation = true;
  bool include_event_invariants = true;

  // --- OccEvent orbit printing options ---
  bool include_elements = false;
  bool include_invariant_group = false;
  bool include_equivalence_map = false;
  xtal::SymInfoOptions sym_info_options;
};

}  // namespace occ_events

// OccPosition

jsonParser &to_json(
    occ_events::OccPosition const &pos, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> systemm);

void from_json(occ_events::OccPosition &pos, jsonParser const &json,
               occ_events::OccSystem const &system);

template <>
struct jsonConstructor<occ_events::OccPosition> {
  static occ_events::OccPosition from_json(jsonParser const &json,
                                           occ_events::OccSystem const &system);
};

void parse(InputParser<occ_events::OccPosition> &parser,
           occ_events::OccSystem const &system);

// OccTrajectory

jsonParser &to_json(
    occ_events::OccTrajectory const &traj, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system);

void from_json(occ_events::OccTrajectory &traj, jsonParser const &json,
               occ_events::OccSystem const &system);

template <>
struct jsonConstructor<occ_events::OccTrajectory> {
  static occ_events::OccTrajectory from_json(
      jsonParser const &json, occ_events::OccSystem const &system);
};

void parse(InputParser<occ_events::OccTrajectory> &parser,
           occ_events::OccSystem const &system);

// OccEvent

jsonParser &to_json(
    occ_events::OccEvent const &event, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system,
    occ_events::OccEventOutputOptions const &options =
        occ_events::OccEventOutputOptions());

/// \brief OccEvent orbit printing
jsonParser &to_json(
    std::set<occ_events::OccEvent> const &orbit, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system,
    std::shared_ptr<occ_events::SymGroup const> const &factor_group,
    std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep,
    occ_events::OccEventOutputOptions const &options =
        occ_events::OccEventOutputOptions());

void from_json(occ_events::OccEvent &event, jsonParser const &json,
               occ_events::OccSystem const &system);

template <>
struct jsonConstructor<occ_events::OccEvent> {
  static occ_events::OccEvent from_json(jsonParser const &json,
                                        occ_events::OccSystem const &system);
};

void parse(InputParser<occ_events::OccEvent> &parser,
           occ_events::OccSystem const &system);

// OccEventInvariants

jsonParser &to_json(occ_events::OccEventInvariants const &invariants,
                    jsonParser &json);

}  // namespace CASM

#endif
