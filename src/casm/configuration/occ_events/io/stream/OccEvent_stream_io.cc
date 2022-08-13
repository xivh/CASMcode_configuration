#include "casm/configuration/occ_events/io/stream/OccEvent_stream_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace occ_events {

// defined in OccEventCounterStateInfo_json_io.cc:
namespace OccEventCounterStateInfoPrinter_impl {

void _print_occ(Log &log, clust::IntegralCluster const &cluster,
                std::vector<int> const &occ, OccSystem const &system);

void _print_traj(Log &log, OccPosition const &pos0, OccPosition const &pos1,
                 OccSystem const &system);

}  // namespace OccEventCounterStateInfoPrinter_impl

void OccEventPrinter::operator()(OccEvent const &event) {
  using namespace OccEventCounterStateInfoPrinter_impl;
  auto &log = m_log;
  auto cluster_occupation = make_cluster_occupation(event);
  auto const &cluster = cluster_occupation.first;
  auto const &occ_init = cluster_occupation.second[0];
  auto const &occ_final = cluster_occupation.second[1];

  if (occ_init.size()) {
    log << "occ_init: ";
    _print_occ(log, cluster, occ_init, m_system);
    log << std::endl;
  }
  if (occ_final.size()) {
    log << "occ_final: ";
    _print_occ(log, cluster, occ_final, m_system);
    log << std::endl;
  }
  if (event.size()) {
    log << "trajectory: " << std::endl;
    for (auto const &traj : event) {
      OccPosition const &pos0 = traj.position[0];
      OccPosition const &pos1 = traj.position[1];
      _print_traj(log, pos0, pos1, m_system);
      log << std::endl;
    }
  }
}

}  // namespace occ_events
}  // namespace CASM
