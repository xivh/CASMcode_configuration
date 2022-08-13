#ifndef CASM_occ_events_OccEvent_stream_io
#define CASM_occ_events_OccEvent_stream_io

#include "casm/casm_io/Log.hh"

namespace CASM {
class Log;

namespace occ_events {
class OccEvent;
struct OccSystem;

/// \brief Print OccEventCounterStateInfo to stream
struct OccEventPrinter {
  OccEventPrinter(OccSystem const &system, CASM::Log const &log = CASM::log())
      : m_system(system), m_log(log) {}

  void operator()(OccEvent const &event);

 private:
  OccSystem const &m_system;
  Log m_log;
};

/// \brief Print IntegralCluster to stream, using default
/// Printer<IntegralCluster>
std::ostream &operator<<(
    std::ostream &sout,
    std::pair<OccEvent const &, OccSystem const &> event_and_system) {
  OccEventPrinter p(event_and_system.second, CASM::Log(sout));
  p(event_and_system.first);
  return sout;
}

}  // namespace occ_events
}  // namespace CASM

#endif
