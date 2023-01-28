#include "casm/configuration/occ_events/io/stream/OccEventCounter_stream_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/occ_events/OccEventCounter.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace occ_events {
namespace OccEventCounterStateInfoPrinter_impl {

void _print_occ(Log &log, clust::IntegralCluster const &cluster,
                std::vector<int> const &occ, OccSystem const &system) {
  log << "[";
  for (Index i = 0; i < cluster.size(); ++i) {
    if (i != 0) {
      log << ", ";
    }
    log << occ[i];
  }
  log << "] == [";
  for (Index i = 0; i < cluster.size(); ++i) {
    if (i != 0) {
      log << ", ";
    }
    log << "\"" << system.get_orientation_name(cluster[i], occ[i]) << "\"";
  }
  log << "]";
}

void _print_traj(Log &log, OccPosition const &pos0, OccPosition const &pos1,
                 OccSystem const &system) {
  jsonParser json;
  log << "[" << to_json(pos0.integral_site_coordinate, json) << ", "
      << pos0.occupant_index;
  if (pos0.is_atom) {
    log << ", " << pos0.atom_position_index;
    log << "] ";
  } else {
    log << "]";
  }
  log << " -> ";
  log << "[" << to_json(pos1.integral_site_coordinate, json) << ", "
      << pos1.occupant_index;
  if (pos0.is_atom) {
    log << ", " << pos1.atom_position_index;
    log << "] ";
  } else {
    log << "]";
  }

  log << "(";
  log << system.get_orientation_name(pos0);
  if (pos0.is_atom) {
    log << ".atom[" << pos0.atom_position_index
        << "]=" << system.get_atom_name(pos0);
  }
  log << ")";
  log << " -> ";
  if (pos0 == pos1) {
    log << "(no change)";
  } else {
    log << "(";
    log << system.get_orientation_name(pos1);
    if (pos0.is_atom) {
      log << ".atom[" << pos1.atom_position_index
          << "]=" << system.get_atom_name(pos1);
    }
    log << ")";
  }
}

}  // namespace OccEventCounterStateInfoPrinter_impl

void OccEventCounterStateInfoPrinter::operator()(
    OccEventCounterStateInfo const &info) {
  using namespace OccEventCounterStateInfoPrinter_impl;
  auto &log = m_log;
  log << "---" << std::endl;
  auto const &cluster = info.cluster;
  auto const &occ_init = info.occ_init;
  auto const &occ_final = info.occ_final;

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
  if (info.position_init.size()) {
    log << "trajectory: " << std::endl;
    for (Index i = 0; i < info.position_init.size(); ++i) {
      OccPosition const &pos0 = info.position_init[i];
      OccPosition const &pos1 = info.position_final[i];
      _print_traj(log, pos0, pos1, m_system);
      log << std::endl;
    }
  }
  if (info.fails.empty()) {
    log << ">> allowed" << std::endl;
  } else {
    log << ">> not allowed, due to " << info.fails << std::endl;
  }
  log << "---" << std::endl;
}

}  // namespace occ_events
}  // namespace CASM
