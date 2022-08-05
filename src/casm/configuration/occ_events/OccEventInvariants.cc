#include "casm/configuration/occ_events/OccEventInvariants.hh"

#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace occ_events {

namespace {

// Only uses first position in each trajectory
std::map<std::string, Index> make_molecule_count(OccEvent const &event,
                                                 OccSystem const &system) {
  std::map<std::string, Index> molecule_count;
  for (auto const &traj : event) {
    if (!traj.position.size()) {
      throw std::runtime_error(
          "Error in `make_molecule_count`: Empty trajectory");
    }
    std::string molecule_name = system.get_name(traj.position[0]);
    auto it = molecule_count.find(molecule_name);
    if (it == molecule_count.end()) {
      molecule_count[molecule_name] = 1;
    } else {
      it->second++;
    }
  }
  return molecule_count;
}

// Determines coordinate for each occupant (site coordinate for Mol, atom
// position for component atoms of molecules, skips resevoir positions) and
// keeps those that are not approximately equal
std::vector<Eigen::Vector3d> make_unique_coordinates(OccEvent const &event,
                                                     OccSystem const &system) {
  double tol = system.prim->lattice().tol();

  std::vector<Eigen::Vector3d> unique_coordinates;
  for (auto const &traj : event) {
    for (auto const &pos : traj.position) {
      if (pos.is_in_resevoir) {
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

  // calculate distances between points (excludes resevoir positions)
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
  m_molecule_count = make_molecule_count(event, system);
  m_distances = make_sorted_distances(event, system);
}

/// Number of elements in the OccEvent
int OccEventInvariants::size() const { return m_size; }

/// Distances between sites in the OccEvent
std::vector<double> const &OccEventInvariants::distances() const {
  return m_distances;
}

/// Count of each molecule type involved in the OccEvent
std::map<std::string, Index> OccEventInvariants::molecule_count() const {
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
  if (A.molecule_count() != B.molecule_count()) {
    return A.molecule_count() < B.molecule_count();
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
  return false;
}

}  // namespace occ_events
}  // namespace CASM
