#ifndef CASM_occ_events_OccEventCounter_stream_io
#define CASM_occ_events_OccEventCounter_stream_io

#include "casm/casm_io/Log.hh"

namespace CASM {
namespace occ_events {
struct OccEventCounterStateInfo;
struct OccSystem;

/// \brief Print OccEventCounterStateInfo to stream
struct OccEventCounterStateInfoPrinter {
  OccEventCounterStateInfoPrinter(OccSystem const &system,
                                  CASM::Log const &log = CASM::log())
      : m_system(system), m_log(log) {}

  void operator()(OccEventCounterStateInfo const &info);

 private:
  OccSystem const &m_system;
  Log m_log;
};

}  // namespace occ_events
}  // namespace CASM

#endif
