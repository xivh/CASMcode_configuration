#ifndef CASM_occ_events_OccEventInvariants
#define CASM_occ_events_OccEventInvariants

#include <map>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {
class BasicStructure;
}

namespace occ_events {

class OccEvent;
struct OccSystem;

class OccEventInvariants {
 public:
  /// \brief Construct and calculate OccEvent invariants
  explicit OccEventInvariants(OccEvent const &event, OccSystem const &system);

  /// Number of elements in the OccEvent
  int size() const;

  /// Sorted distances between sites in the OccEvent
  /// - Excludes sites with sublattice == -1 (which is used to indicate
  ///   the resevoir)
  std::vector<double> const &distances() const;

  /// Count of each molecule type involved in the OccEvent
  std::map<std::string, Index> molecule_count() const;

 private:
  /// Number of elements in the OccEvent
  int m_size;

  /// Distances between sites in the OccEvent
  std::vector<double> m_distances;

  /// Count of each molecule type involved in the OccEvent
  std::map<std::string, Index> m_molecule_count;
};

/// Check if OccEventInvariants are equal
bool almost_equal(OccEventInvariants const &A, OccEventInvariants const &B,
                  double tol);

/// Compare OccEventInvariants
bool compare(OccEventInvariants const &A, OccEventInvariants const &B,
             double tol);

}  // namespace occ_events
}  // namespace CASM

#endif
