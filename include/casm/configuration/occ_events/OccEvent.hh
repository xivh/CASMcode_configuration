#ifndef CASM_occ_events_OccEvent
#define CASM_occ_events_OccEvent

#include <vector>

#include "casm/configuration/occ_events/OccTrajectory.hh"
#include "casm/configuration/occ_events/definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace occ_events {

/// OccEvent describes changes in occupation
///
/// - It can be used to find orbits of symmetrically equivalent occupation
///   events
/// - It can be used to find the invariant group of occupation changes for
///   generating local cluster expansions
class OccEvent : public Comparisons<xtal::Translatable<CRTPBase<OccEvent>>> {
 public:
  typedef std::vector<OccTrajectory>::iterator iterator;
  typedef std::vector<OccTrajectory>::const_iterator const_iterator;
  typedef std::vector<OccTrajectory>::size_type size_type;

  explicit OccEvent();

  explicit OccEvent(std::initializer_list<OccTrajectory> elements);

  explicit OccEvent(std::vector<OccTrajectory> const &elements);

  explicit OccEvent(std::vector<OccTrajectory> &&elements);

  template <typename Iterator>
  OccEvent(Iterator begin, Iterator end);

  /// Access vector of elements
  std::vector<OccTrajectory> &elements();

  /// const Access vector of elements
  const std::vector<OccTrajectory> &elements() const;

  /// Iterator to first element in the OccEvent
  iterator begin();

  /// Iterator to first element in the OccEvent
  const_iterator begin() const;

  /// Iterator to the past-the-last element in the OccEvent
  iterator end();

  /// Iterator to the past-the-last element in the OccEvent
  const_iterator end() const;

  /// Iterator to first element in the OccEvent
  const_iterator cbegin() const;

  /// Iterator to the past-the-last element in the OccEvent
  const_iterator cend() const;

  /// Number of elements in the OccEvent
  size_type size() const;

  /// Access an element in the OccEvent by index
  OccTrajectory &operator[](size_type index);

  /// Access an element in the OccEvent by index
  OccTrajectory const &operator[](size_type index) const;

  /// Access an element in the OccEvent by index
  OccTrajectory &element(size_type index);

  /// Access a UnitCellCoord in the OccEvent by index
  OccTrajectory const &element(size_type index) const;

  /// Translate the OccEvent by a UnitCell translation
  OccEvent &operator+=(xtal::UnitCell trans);

  /// Basic comparison (does not consider reverse)
  bool operator<(OccEvent const &B) const;

 private:
  std::vector<OccTrajectory> m_element;
};

OccEvent &sort(OccEvent &event);

OccEvent copy_sort(OccEvent event);

OccEvent &reverse(OccEvent &event);

OccEvent copy_reverse(OccEvent event);

/// \brief Put event into standardized form with regard to permutation/reversal
OccEvent &standardize(OccEvent &occ_event);

/// \brief Make an OccEvent from before and after positions
OccEvent make_occevent(std::vector<OccPosition> const &position_before,
                       std::vector<OccPosition> const &position_after);

/// \brief Make clust::IntegralCluster from first position in each trajectory
clust::IntegralCluster make_cluster(OccEvent const &event);

/// \brief Make {cluster, {occ_init, occ_final}} from OccEvent
std::pair<clust::IntegralCluster, std::vector<std::vector<int>>>
make_cluster_occupation(OccEvent const &event);

/// \brief Apply SymOp to OccEvent
OccEvent &apply(OccEventRep const &rep, OccEvent &occ_event);

/// \brief Apply SymOp to OccEvent
OccEvent copy_apply(OccEventRep const &rep, OccEvent occ_event);

// --- Implementation ---

template <typename Iterator>
OccEvent::OccEvent(Iterator begin, Iterator end) : m_element(begin, end) {}

}  // namespace occ_events
}  // namespace CASM

#endif
