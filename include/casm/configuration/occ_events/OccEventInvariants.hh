#ifndef CASM_occ_events_OccEventInvariants
#define CASM_occ_events_OccEventInvariants

#include <map>
#include <set>
#include <vector>

#include "casm/configuration/occ_events/misc/LexicographicalCompare.hh"
#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

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
  ///   the reservoir)
  std::vector<double> const &distances() const;

  /// Count of each molecule type involved in the OccEvent,
  ///     in initial / final states if different
  std::set<Eigen::VectorXi, LexicographicalCompare> const &molecule_count()
      const;

 private:
  /// Number of elements in the OccEvent
  int m_size;

  /// Distances between sites in the OccEvent
  std::vector<double> m_distances;

  /// Count of each molecule type involved in the OccEvent,
  ///     in initial / final states if different
  std::set<Eigen::VectorXi, LexicographicalCompare> m_molecule_count;
};

/// Check if OccEventInvariants are equal
bool almost_equal(OccEventInvariants const &A, OccEventInvariants const &B,
                  double tol);

/// Compare OccEventInvariants
bool compare(OccEventInvariants const &A, OccEventInvariants const &B,
             double tol);

/// \brief Compare std::pair<OccEventInvariants, OccEvent>
struct CompareOccEvent_f {
  typedef std::pair<OccEventInvariants, OccEvent> pair_type;

  CompareOccEvent_f(double _xtal_tol) : xtal_tol(_xtal_tol) {}

  bool operator()(pair_type const &A, pair_type const &B) const;

  double xtal_tol;
};

}  // namespace occ_events
}  // namespace CASM

#endif
