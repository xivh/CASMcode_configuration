#include "casm/configuration/occ_events/OccEvent.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/OccTrajectory.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace occ_events {

OccEvent::OccEvent() {}

OccEvent::OccEvent(std::initializer_list<OccTrajectory> elements)
    : m_element(elements) {}

OccEvent::OccEvent(std::vector<OccTrajectory> const &elements)
    : m_element(elements) {}

OccEvent::OccEvent(std::vector<OccTrajectory> &&elements)
    : m_element(std::move(elements)) {}

/// \brief Access vector of elements
std::vector<OccTrajectory> &OccEvent::elements() { return m_element; }

/// \brief const Access vector of elements
std::vector<OccTrajectory> const &OccEvent::elements() const {
  return m_element;
}

/// \brief Iterator to first element in the OccEvent
typename OccEvent::iterator OccEvent::begin() { return elements().begin(); }

/// \brief Iterator to first element in the OccEvent
typename OccEvent::const_iterator OccEvent::begin() const {
  return elements().begin();
}

/// \brief Iterator to the past-the-last element in the OccEvent
typename OccEvent::iterator OccEvent::end() { return elements().end(); }

/// \brief Iterator to the past-the-last element in the OccEvent
typename OccEvent::const_iterator OccEvent::end() const {
  return elements().end();
}

/// \brief Iterator to first element in the OccEvent
typename OccEvent::const_iterator OccEvent::cbegin() const {
  return elements().cbegin();
}

/// \brief Iterator to the past-the-last element in the OccEvent
typename OccEvent::const_iterator OccEvent::cend() const {
  return elements().cend();
}

/// \brief Number of elements in the OccEvent
typename OccEvent::size_type OccEvent::size() const {
  return elements().size();
}

/// \brief Access an element in the OccEvent by index
OccTrajectory &OccEvent::operator[](size_type index) { return element(index); }

/// \brief Access an element in the OccEvent by index
OccTrajectory const &OccEvent::operator[](size_type index) const {
  return element(index);
}

/// \brief Access an element in the OccEvent by index
OccTrajectory &OccEvent::element(size_type index) { return elements()[index]; }

/// \brief Access a UnitCellCoord in the OccEvent by index
OccTrajectory const &OccEvent::element(size_type index) const {
  return elements()[index];
}

/// \brief Translate the OccEvent by a UnitCell translation
OccEvent &OccEvent::operator+=(xtal::UnitCell trans) {
  for (auto it = this->begin(); it != this->end(); ++it) {
    *it += trans;
  }
  return *this;
}

/// Basic comparison (does not consider reverse)
bool OccEvent::operator<(OccEvent const &B) const {
  if (this->size() != B.size()) {
    return this->size() < B.size();
  }
  return this->elements() < B.elements();
}

OccEvent &sort(OccEvent &event) {
  std::sort(event.elements().begin(), event.elements().end());
  return event;
}

OccEvent copy_sort(OccEvent event) { return sort(event); }

OccEvent &reverse(OccEvent &event) {
  for (auto &e : event.elements()) {
    std::reverse(e.position.begin(), e.position.end());
  }
  return event;
}

OccEvent copy_reverse(OccEvent event) { return reverse(event); }

/// \brief Put event into standardized form with regard to permutation/reversal
OccEvent &standardize(OccEvent &occ_event) {
  OccEvent reverse_occ_event = copy_reverse(occ_event);
  sort(occ_event);
  sort(reverse_occ_event);
  if (reverse_occ_event < occ_event) {
    occ_event = reverse_occ_event;
  }
  return occ_event;
}

/// \brief Make an OccEvent from before and after positions
///
/// \param position_before The initial positions
/// \param position_final The final positions
///
/// Note:
/// - This makes an OccEvent consiting of
///   OccTrajectory({position_before[i], position_after[i]})
///   for each i.
/// - This does no validity checking
OccEvent make_occevent(std::vector<OccPosition> const &position_before,
                       std::vector<OccPosition> const &position_after) {
  std::vector<OccTrajectory> traj(position_before.size());
  for (Index i = 0; i < position_before.size(); ++i) {
    traj[i].position.emplace_back(position_before[i]);
    traj[i].position.emplace_back(position_after[i]);
  }
  return OccEvent(std::move(traj));
}

/// \brief Make IntegralCluster from first position in each trajectory
///
/// - Skip any position that `is_in_reservoir`
/// - Throw for any OccTrajectory with empty position vector
/// - Result is sorted
clust::IntegralCluster make_cluster(OccEvent const &event) {
  std::set<xtal::UnitCellCoord> integral_sites;
  for (OccTrajectory const &traj : event) {
    if (!traj.position.size()) {
      throw std::runtime_error(
          "Error in `make_cluster(OccEvent const&)`: Empty trajectory");
    }
    OccPosition const &pos = traj.position[0];
    if (pos.is_in_reservoir) {
      continue;
    }
    integral_sites.insert(pos.integral_site_coordinate);
  }
  clust::IntegralCluster cluster(integral_sites.begin(), integral_sites.end());
  cluster.sort();
  return cluster;
}

/// \brief Make {cluster, {occ_init, occ_final}} from OccEvent
///
/// Notes:
/// - Throws if any OccTrajectory is not size 2
/// - Throws if any inconsistency in site occupation
std::pair<clust::IntegralCluster, std::vector<std::vector<int>>>
make_cluster_occupation(OccEvent const &event) {
  std::pair<clust::IntegralCluster, std::vector<std::vector<int>>> result(
      clust::IntegralCluster(),
      std::vector<std::vector<int>>(2, std::vector<int>()));
  if (!event.size()) {
    return result;
  }
  for (auto const &traj : event) {
    if (traj.position.size() != 2) {
      throw std::runtime_error(
          "Error in make_cluster_occupation: a trajectory does not have 2 "
          "positions");
    }
  }
  result.first = make_cluster(event);
  auto const &cluster = result.first;
  result.second =
      std::vector<std::vector<int>>(2, std::vector<int>(cluster.size(), -1));
  auto &cluster_occupation = result.second;

  for (OccTrajectory const &traj : event) {
    for (Index s = 0; s < 2; ++s) {
      OccPosition const &pos = traj.position[s];
      if (pos.is_in_reservoir) {
        continue;
      }
      auto const &site = pos.integral_site_coordinate;
      Index i = find_index(cluster, site);
      if (cluster_occupation[s][i] == -1) {
        cluster_occupation[s][i] = pos.occupant_index;
      } else if (cluster_occupation[s][i] != pos.occupant_index) {
        throw std::runtime_error(
            "Error in make_cluster_occupation: Inconsistent occupancies");
      }
    }
  }
  return result;
}

/// \brief Apply SymOp to OccEvent
OccEvent &apply(OccEventRep const &rep, OccEvent &occ_event) {
  for (auto &traj : occ_event) {
    apply(rep, traj);
  }
  return occ_event;
}

/// \brief Apply SymOp to OccEvent
OccEvent copy_apply(OccEventRep const &rep, OccEvent occ_event) {
  apply(rep, occ_event);
  return occ_event;
}

}  // namespace occ_events
}  // namespace CASM
