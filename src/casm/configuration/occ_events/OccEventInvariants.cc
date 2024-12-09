#include "casm/configuration/occ_events/OccEventInvariants.hh"

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/CASM_math.hh"

// debug
#include "casm/casm_io/container/stream_io.hh"

namespace CASM {
namespace occ_events {

namespace {

// Determines coordinate for each occupant (site coordinate for Mol, atom
// position for component atoms of molecules, skips reservoir positions) and
// keeps those that are not approximately equal
std::vector<Eigen::Vector3d> make_unique_coordinates(OccEvent const &event,
                                                     OccSystem const &system) {
  double tol = system.prim->lattice().tol();

  std::vector<Eigen::Vector3d> unique_coordinates;
  for (auto const &traj : event) {
    for (auto const &pos : traj.position) {
      if (pos.is_in_reservoir) {
        continue;
      }
      Eigen::Vector3d test_coord = system.get_cartesian_coordinate(pos);
      auto begin = unique_coordinates.begin();
      auto end = unique_coordinates.end();
      auto almost_equal_f = [&](Eigen::Vector3d const &existing_coord) {
        return almost_equal(test_coord, existing_coord, tol);
      };
      if (std::find_if(begin, end, almost_equal_f) == end) {
        unique_coordinates.push_back(test_coord);
      }
    }
  }
  return unique_coordinates;
}

// Calculate distance between each unique pair of coordinates, then sorts
std::vector<double> make_sorted_distances(OccEvent const &event,
                                          OccSystem const &system) {
  std::vector<Eigen::Vector3d> unique_coordinates =
      make_unique_coordinates(event, system);

  // calculate distances between points (excludes reservoir positions)
  std::vector<double> distances;
  for (int i = 0; i < unique_coordinates.size(); i++) {
    for (int j = i + 1; j < unique_coordinates.size(); j++) {
      distances.push_back(
          (unique_coordinates[i] - unique_coordinates[j]).norm());
    }
  }
  // and sort
  std::sort(distances.begin(), distances.end());
  return distances;
}

}  // namespace

/// \brief Construct and calculate OccEvent invariants
OccEventInvariants::OccEventInvariants(OccEvent const &event,
                                       OccSystem const &system) {
  m_size = event.size();
  m_distances = make_sorted_distances(event, system);

  // {cluster, {occ_init, occ_final}}
  auto pair = make_cluster_occupation(event);
  Eigen::VectorXi count;
  system.molecule_count(count, pair.first, pair.second[0]);
  m_molecule_count.insert(count);
  system.molecule_count(count, pair.first, pair.second[1]);
  m_molecule_count.insert(count);
}

/// Number of elements in the OccEvent
int OccEventInvariants::size() const { return m_size; }

/// Distances between sites in the OccEvent
std::vector<double> const &OccEventInvariants::distances() const {
  return m_distances;
}

/// Count of each molecule type involved in the OccEvent
std::set<Eigen::VectorXi, LexicographicalCompare> const &
OccEventInvariants::molecule_count() const {
  return m_molecule_count;
}

/// Check if OccEventInvariants are equal
bool almost_equal(OccEventInvariants const &A, OccEventInvariants const &B,
                  double tol) {
  return A.size() == B.size() && A.molecule_count() == B.molecule_count() &&
         std::equal(A.distances().cbegin(), A.distances().cend(),
                    B.distances().cbegin(), [&](double a, double b) {
                      return CASM::almost_equal(a, b, tol);
                    });
}

/// Compare OccEventInvariants
bool compare(OccEventInvariants const &A, OccEventInvariants const &B,
             double tol) {
  // first sort by number of sites in cluster
  if (A.size() != B.size()) {
    return A.size() < B.size();
  }
  // all distances
  for (int i = A.distances().size() - 1; i >= 0; i--) {
    if (CASM::almost_equal(A.distances()[i], B.distances()[i], tol)) {
      continue;
    }
    if (A.distances()[i] < B.distances()[i]) {
      return true;
    }
    if (A.distances()[i] > B.distances()[i]) {
      return false;
    }
  }
  if (A.molecule_count() != B.molecule_count()) {
    // return A.molecule_count() < B.molecule_count();
    return !std::lexicographical_compare(
        A.molecule_count().begin(), A.molecule_count().end(),
        B.molecule_count().begin(), B.molecule_count().end(),
        LexicographicalCompare());
  }
  return false;
}

bool CompareOccEvent_f::operator()(pair_type const &A,
                                   pair_type const &B) const {
  if (compare(A.first, B.first, xtal_tol)) {
    return true;
  }
  if (compare(B.first, A.first, xtal_tol)) {
    return false;
  }
  return A.second < B.second;
}

}  // namespace occ_events
}  // namespace CASM
