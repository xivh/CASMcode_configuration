#ifndef CASM_occ_events_OccEventCounter
#define CASM_occ_events_OccEventCounter

#include <optional>
#include <set>
#include <vector>

#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/occ_counter.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/misc/MultiStepMethod.hh"
#include "casm/container/Counter.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace occ_events {

/// \brief OccEventCounter parameters
///
/// Notes:
/// - These parameters can be used to filter which events are
///   generated.
struct OccEventCounterParameters {
  // --- Cluster checks ---

  std::optional<int> min_cluster_size;
  std::optional<int> max_cluster_size;
  std::optional<int> required_cluster_size;
  std::optional<std::set<int>> excluded_sublattices;
  std::optional<std::set<int>> required_sublattices;

  // --- Initial occupation checks ---

  // Filter events by how many atom/molecule/orientations
  // exist before the event occurs. Note that the counts
  // apply to the event as generated, not the canonical
  // form.

  // Order determined by OccSystem::atom_name_list
  std::optional<Eigen::VectorXi> min_init_atom_count;
  std::optional<Eigen::VectorXi> max_init_atom_count;
  std::optional<Eigen::VectorXi> required_init_atom_count;

  // Order determined by OccSystem::molecule_name_list
  std::optional<Eigen::VectorXi> min_init_molecule_count;
  std::optional<Eigen::VectorXi> max_init_molecule_count;
  std::optional<Eigen::VectorXi> required_init_molecule_count;

  // Order determined by OccSystem::orientation_name_list
  std::optional<Eigen::VectorXi> min_init_orientation_count;
  std::optional<Eigen::VectorXi> max_init_orientation_count;
  std::optional<Eigen::VectorXi> required_init_orientation_count;

  // --- Initial + Final occupation checks ---

  // Order determined by OccSystem::atom_name_list
  std::optional<Eigen::VectorXi> min_final_atom_count;
  std::optional<Eigen::VectorXi> max_final_atom_count;
  std::optional<Eigen::VectorXi> required_final_atom_count;

  // Order determined by OccSystem::molecule_name_list
  std::optional<Eigen::VectorXi> min_final_molecule_count;
  std::optional<Eigen::VectorXi> max_final_molecule_count;
  std::optional<Eigen::VectorXi> required_final_molecule_count;

  // Order determined by OccSystem::orientation_name_list
  std::optional<Eigen::VectorXi> min_final_orientation_count;
  std::optional<Eigen::VectorXi> max_final_orientation_count;
  std::optional<Eigen::VectorXi> required_final_orientation_count;

  // --- Trajectory checks ---

  /// \brief Skip or allow events which can be generated
  ///     using a subcluster because (i) a site is vacant
  ///     before and after the event or (ii) any molecule
  ///     does not change sites, break up, or re-orient.
  bool allow_subcluster_events = false;

  /// \brief Skip or allow events in which molecules break up
  ///     (i.e. do not allow dumbbell interstitials to separate)
  bool do_not_allow_breakup = false;

  // --- debugging ---

  /// \brief If true, allow occ_final < occ_init, otherwise
  ///     this is not allowed to skip generating an event
  ///     and its reverse event
  bool allow_reverse_occ = false;

  /// \brief Current implementation requires this is true.
  ///     Future implementation will allow moves to/from resevoir
  ///     by setting this to false.
  bool require_atom_conservation = true;

  /// \brief If not allowing molecule breakup or moves to/from resevoir,
  ///     requiring molecule conservation can allow skipping some
  ///     unallowed final occupations.
  bool require_molecule_conservation = false;

  /// \brief Skip trajectories in which the atom/molecule type
  ///     changes (should always be true except for debugging
  ///     purposes)
  bool require_chemical_type_conserving_trajectories = true;

  /// \brief Skip or allow events in which indivisible molecules
  ///     break up (should always be true except for debugging
  ///     purposes)
  bool do_not_allow_indivisible_molecule_breakup = true;

  /// \brief Skip events in which a molecule in the resevoir
  ///     remains in the resevoir (for debugging purposes...
  ///     this situation should not occur)
  bool require_no_molecules_remain_in_resevoir = true;
};

/// \brief Data structure used internally by OccEventCounter
struct OccEventCounterData {
  /// \brief Defines OccPosition indices, helps check OccEvents
  std::shared_ptr<OccSystem> system;

  /// \brief All clusters on which OccEvents will be generated
  std::vector<clust::IntegralCluster> prototypes;

  /// \brief Parameters controlling the OccEventCounter
  OccEventCounterParameters params;

  /// \brief Counts over prototype vector. Incremented in
  ///     the outer-most step (step index 3).
  Index prototype_index;

  /// \brief Current cluster on which OccEvents are being generated,
  ///     a copy of prototype[prototype_index].
  clust::IntegralCluster cluster;

  /// \brief Counter generates initial occupation states on
  ///     the cluster. Incremented in second-outer-most step
  ///     (step index 2).
  Counter<std::vector<int>> occ_init_counter;

  /// \brief Counter generates final occupation states on
  ///     the cluster. Incremented in third-outer-most step
  ///     (step index 1).
  Counter<std::vector<int>> occ_final_counter;

  /// \brief OccPosition determined directly from occ_init_counter.
  ///
  /// Notes:
  /// - May be atom positions or molecule positions, depending on
  ///   `params.enumerate_individual_atom_trajectories`
  std::vector<OccPosition> position_init;

  /// \brief A permutation of `position_init`. The trajectories of
  ///     atoms or molecules are generated as position_init[i] ->
  ///     position_init[j]. Permuted in the inner-most step (step
  ///     index 0). Some permutations are allowable and some are
  ///     not because they result in atom/molecule type change,
  ///     which is checked during the method. A parameter
  ///     (`allow_subcluster_events`) determines whether permutations
  ///     which leave some sites unchanged (subcluster events) are
  ///     allowed.
  ///
  std::vector<OccPosition> position_final;

  /// \brief Current OccEvent
  OccEvent occ_event;
};

/// \brief Count over potential OccEvent
///
/// Notes:
/// - Given a vector of clust::IntegralCluster, this method iterates over
///   each cluster and generates all possible OccEvent (under user-specified
///   conditions).
/// - The method consists of 4 steps in a nested loop:
///   - iterate over cluster prototypes (outer-most step)
///   - iterate over initial cluster occupation
///   - iterate over final cluster occupation
///   - iterate over initial->final trajectory permutations (inner-most step)
/// - At each step a number of criteria, specified by OccEventCounterParameters,
///   are checked to skip invalid or undesired OccEvent.
class OccEventCounter {
 public:
  /// \brief Constructor
  ///
  /// \param system, OccSystem used to define and check OccEvent
  /// \param clusters, Vector of underlying cluster orbit prototypes
  ///     on which OccEvent should be generated.
  /// \param params, Options controlling the events generated.
  ///
  OccEventCounter(std::shared_ptr<OccSystem> const &system,
                  std::vector<clust::IntegralCluster> const &prototypes,
                  OccEventCounterParameters const &params);

  /// \brief Advance to the next allowed OccEvent
  bool advance();

  /// \brief Current OccEvent
  OccEvent const &value() const;

  /// \brief True if counting is finished
  bool is_finished() const;

 private:
  std::shared_ptr<OccEventCounterData> m_data;
  std::unique_ptr<MultiStepMethod<OccEventCounterData>> m_stepper;
};

}  // namespace occ_events
}  // namespace CASM

#endif
