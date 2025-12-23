#ifndef CASM_group_subgroups
#define CASM_group_subgroups

#include <algorithm>
#include <map>

#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/definitions.hh"

namespace CASM {
namespace group {

typedef std::set<Index> SubgroupIndices;
typedef std::set<SubgroupIndices> SubgroupOrbit;

/// \brief Return all cyclic subgroups
template <typename ElementType>
std::set<SubgroupOrbit> make_cyclic_subgroups(Group<ElementType> const &group);

/// \brief Return all subgroups
template <typename ElementType>
std::set<SubgroupOrbit> make_all_subgroups(Group<ElementType> const &group);

/// \brief Make the invariant subgroup for each orbit element, as
///     indices of group elements
template <typename GroupElementType>
std::vector<SubgroupIndices> make_invariant_subgroups(
    std::vector<std::vector<Index>> const &equivalence_map,
    Group<GroupElementType> const &group);

/// \brief Utility class generating a subgroup from known generators
template <typename ElementType>
struct MakeSubgroupFromGenerators {
  /// \brief Elements generated upon closure
  std::set<Index> indices;

  /// \brief True if element in head group is member of subgroup
  std::vector<bool> member;

  /// \brief The generating elements
  std::set<Index> generators;

  /// \brief Constructor
  ///
  /// \param group The head group, with multiplication table to use
  /// \param generators Indices of known subgroup generators
  MakeSubgroupFromGenerators(Group<ElementType> const &group)
      : member(group.element.size(), false) {}

  /// \brief Constructor
  ///
  /// \param group The head group, with multiplication table to use
  /// \param generators Indices of known subgroup generators
  MakeSubgroupFromGenerators(Group<ElementType> const &group,
                             std::set<Index> const &_generators)
      : member(group.element.size(), false), generators(_generators) {
    /// Add identity to the result:
    indices.insert(0);
    member[0] = true;

    /// Queue of products - start with identity
    std::set<Index> queue;
    queue.insert(0);

    _close(group, queue);
  }

  void add_generator(Group<ElementType> const &group, Index k) {
    if (member[k]) {
      return;
    }
    generators.insert(k);

    std::set<Index> queue;

    if (indices.size()) {
      /// Products of the existing elements and the new generator
      /// to get new elements which are added to the queue
      Index prod;
      for (Index i : indices) {
        prod = group.mult(i, k);
        if (!indices.count(prod)) {
          queue.insert(prod);
        }
      }

    } else {
      /// Add identity to the result:
      indices.insert(0);
      member[0] = true;

      /// Queue of products - start with identity
      queue.insert(0);
    }

    _close(group, queue);
  }

 private:
  void _close(Group<ElementType> const &group, std::set<Index> &queue) {
    /// Products of the new elements and all generators
    Index prod;
    do {
      auto begin = queue.begin();
      Index x = *begin;
      queue.erase(*begin);

      // For each generator:
      for (auto j : generators) {
        // Check x * j
        prod = group.mult(x, j);
        if (!member[prod]) {
          indices.insert(prod);
          member[prod] = true;
          queue.insert(prod);
        }
        // Check x * j^-1
        prod = group.mult(x, group.inv(j));
        if (!member[prod]) {
          indices.insert(prod);
          member[prod] = true;
          queue.insert(prod);
        }
      }

    } while (queue.size());
  }
};

/// \brief Return true if `a` is a subset of `b` (may be equal)
inline bool is_subset_of(std::set<Index> const &a, std::set<Index> const &b) {
  if (a.size() > b.size()) {
    return false;
  }

  // Both sets are ordered; use std::includes to check if `this` contains
  // all elements of `other`.
  return std::includes(b.begin(), b.end(), a.begin(), a.end());
}

/// \brief Return true if `a` is a proper subset of `b` (may not be equal)
inline bool is_proper_subset_of(std::set<Index> const &a,
                                std::set<Index> const &b) {
  if (a.size() == b.size()) {
    return false;
  }
  return is_subset_of(a, b);
}

/// \brief Utility class finding a minimal set of generators
template <typename ElementType>
struct MakeMaximalCyclicSubgroups {
  /// Cyclic subgroups which are not contained by another cyclic subgroup
  std::vector<std::set<Index>> maximal_cyclic_subgroups;

  /// Elements which generate the corresponding maximal cyclic subgroup
  std::vector<Index> generators;

  /// \brief Constructor
  ///
  /// \param group The head group, with multiplication table to use
  /// \param indices Indices of subset elements
  MakeMaximalCyclicSubgroups(Group<ElementType> const &group,
                             std::set<Index> const &indices) {
    for (Index x : indices) {
      std::set<Index> new_cycle;
      Index prod = x;
      while (new_cycle.insert(prod).second) {
        prod = group.mult(x, prod);
      }

      _update(x, new_cycle);
    }
  }

 private:
  void _update(Index new_generator, std::set<Index> const &new_cycle) {
    Index i = 0;
    while (i != maximal_cyclic_subgroups.size()) {
      auto const &existing = maximal_cyclic_subgroups[i];
      if (is_subset_of(new_cycle, existing)) {
        return;
      }
      if (is_proper_subset_of(existing, new_cycle)) {
        if (i != maximal_cyclic_subgroups.size() - 1) {
          generators[i] = std::move(generators.back());
          maximal_cyclic_subgroups[i] = maximal_cyclic_subgroups.back();
        }
        generators.pop_back();
        maximal_cyclic_subgroups.pop_back();
      } else {
        ++i;
      }
    }
    generators.push_back(new_generator);
    maximal_cyclic_subgroups.push_back(new_cycle);
  }
};

/// \brief Utility class finding a minimal set of generators
template <typename ElementType>
struct MakeMinimalSubsetGenerators {
  /// \brief A minimal set of generators of the subset
  std::set<Index> generators;

  /// Elements generated upon closure - may include elements not in the
  /// original subset if the subset is not a group
  std::set<Index> indices;

  /// \brief Constructor
  ///
  /// \param group The head group, with multiplication table to use
  /// \param indices Indices of subset elements
  MakeMinimalSubsetGenerators(Group<ElementType> const &group,
                              std::set<Index> const &subset_indices) {
    MakeMaximalCyclicSubgroups f(group, subset_indices);

    // Make a lookup table for conjugacy classes included in the new subgroup
    std::vector<bool> classes_included;
    for (Index i : subset_indices) {
      Index cc = group.class_of(i);
      if (cc >= classes_included.size()) {
        classes_included.resize(cc + 1);
      }
      classes_included[cc] = false;
    }

    // Loop over maximal cyclic subgroup generators,
    // adding the ones that add new conjugacy classes first
    MakeSubgroupFromGenerators x(group);
    for (Index i = 0; i < f.generators.size(); ++i) {
      // Update classes included in the subgroup already
      for (Index j : x.indices) {
        classes_included[group.class_of(j)] = true;
      }

      // Select best next generator
      std::optional<Index> best_generator = std::nullopt;
      Index best_n_new_class = 0;
      for (Index k_index = 0; k_index < f.generators.size(); ++k_index) {
        Index k = f.generators[k_index];
        if (x.generators.count(k)) {
          continue;
        }
        if (x.indices.count(k)) {
          continue;
        }
        Index n_new_class = 0;
        for (Index j : f.maximal_cyclic_subgroups[k_index]) {
          if (!classes_included[group.class_of(j)]) {
            n_new_class++;
          }
        }
        if (!best_generator.has_value()) {
          best_generator = k;
          best_n_new_class = n_new_class;
        } else if (n_new_class > best_n_new_class) {
          best_generator = k;
          best_n_new_class = n_new_class;
        }
      }

      if (!best_generator.has_value()) {
        break;
        // throw std::runtime_error(
        //     "Error in MakeMinimalSubsetGenerators: logic error");
      }

      x.add_generator(group, *best_generator);
    }
    generators = x.generators;
    indices = x.indices;
  }
};

template <typename ElementType>
struct MakeAllSubgroupsFromGenerators {
  /// \brief Map type: (indices) -> (generators)
  typedef std::map<std::set<Index>, std::set<Index>> map_type;

  /// \brief Map iterator type
  typedef map_type::iterator map_iterator;

  /// \brief All subgroups, as the map (indices) -> (generators)
  map_type subgroups;

  /// \brief Constructor
  ///
  /// This assumes this subset is a group.
  MakeAllSubgroupsFromGenerators(Group<ElementType> const &group,
                                 std::set<Index> const &indices) {
    // std::cout << "begin MakeAllSubgroupsFromGenerators" << std::endl;

    // Results: (indices -> generators)
    Index group_size = group.element.size();

    // Perform a depth first search for combinations of elements to use as
    // subgroup generators

    // Example:
    // - brackets indicate the set of generators for a subgroup
    // - don't need to check {identity}
    //
    // {1} -> {1,2} -> {1,2,3}  (stop if full group generated)
    //              -> {1,2,4} -> {1,2,4,5} (stop...)
    //              -> {1,2,5} (stop...)
    //     -> {1,3} (stop...)
    //     -> {1,4} (stop...)
    //     -> {1,5} (stop...)
    //     -> {1,6} (stop...)
    // {2} -> {2,3} (stop...)
    // {3} -> {3,4} -> ...
    //

    // typedefs
    typedef MakeSubgroupFromGenerators<ElementType> node_type;
    typedef std::set<Index>::const_iterator member_iterator;

    // Get the first element in the subset
    member_iterator member_it = indices.begin();
    if (*member_it != 0) {
      throw std::runtime_error(
          "Error in Subset::all_subgroups: No identity element");
    }
    if (group_size == 1) {
      subgroups.emplace(std::set<Index>{0}, std::set<Index>{0});
      // std::cout << "# subgroups found: " << subgroups.size() << std::endl;
      // std::cout << "end MakeAllSubgroupsFromGenerators" << std::endl;
      return;
    }

    // This is a vector of nodes on the depth-first search
    std::vector<node_type> nodes;

    // This is a vector of iterators into the members of the subset so that
    // we know the next node to search at each level
    std::vector<member_iterator> iters;

    ++member_it;
    iters.push_back(member_it);

    auto print_position = [&]() {
      std::cout << "---" << std::endl;
      std::cout << "nodes: ";
      for (auto const &n : nodes) {
        std::cout << "{";
        for (auto const &i : n.generators) {
          std::cout << i << ",";
        }
        std::cout << "} ";
      }
      std::cout << std::endl;
      std::cout << "iters: ";
      for (auto const &it : iters) {
        if (it != indices.end()) {
          std::cout << *it << " ";
        } else {
          std::cout << "end ";
        }
      }
      std::cout << std::endl << std::endl;
    };

    // Begin
    // std::cout << "begin loop" << std::endl;
    std::pair<map_iterator, bool> emplace_res;
    while (true) {
      // Backtrack over exhausted iterators before preparing the next node.
      // If the last iterator is `indices.end()` we should pop it (and the
      // corresponding node if present) and continue backtracking. This
      // prevents dereferencing `indices.end()` below when the subset has
      // only the identity element (or when deeper levels are exhausted).
      do {
        // print_position();
        while (nodes.size() > 0 && iters.size() > 0 &&
               iters.back() == indices.end()) {
          // std::cout << "Backtracking over exhausted iterator" << std::endl;
          nodes.pop_back();
          iters.pop_back();
          ++(iters.back());
          // print_position();
        }
        if (nodes.size() == 0) {
          break;
        }

        auto const &last_node = nodes.back();
        auto next_it = iters.back();
        // std::cout << "Checking if generator " << *next_it
        //           << " is already in subgroup indices" << std::endl;
        if (!last_node.indices.count(*next_it)) {
          // std::cout << "Generator not already in subgroup indices,
          // proceeding"
          //           << std::endl;
          break;
        }

        // std::cout << "Generator already in subgroup indices, advancing
        // iterator"
        //           << std::endl;
        ++(iters.back());
      } while (true);

      // Check for termination:
      if (nodes.size() == 0 && iters.back() == indices.end()) {
        // std::cout << "All iterators exhausted, ending search" << std::endl;
        // std::cout << std::endl;
        break;
      }

      // Create next node
      if (nodes.size() == 0) {
        // std::cout << "New root node" << std::endl;
        nodes.emplace_back(group);
      } else {
        // std::cout << "Add to existing node" << std::endl;
        auto const &last_node = nodes.back();
        nodes.emplace_back(last_node);
      }
      auto &next_node = nodes.back();
      // std::cout << "Existing generators: ";
      // for (auto const &g : next_node.generators) {
      //   std::cout << g << " ";
      // }
      // std::cout << std::endl;
      // std::cout << "Indices: ";
      // for (auto const &i : next_node.indices) {
      //   std::cout << i << " ";
      // }
      // std::cout << std::endl;

      // Add the new generator
      auto next_it = iters.back();
      if (nodes.size() == 1) {
        std::cout << "Adding generator " << *next_it << std::endl;
      }
      next_node.add_generator(group, *next_it);

      // std::cout << "Indices: ";
      // for (auto const &i : next_node.indices) {
      //   std::cout << i << " ";
      // }
      // std::cout << std::endl;

      // Try adding subgroup
      emplace_res = subgroups.emplace(next_node.indices, next_node.generators);

      // if (emplace_res.second) {
      //   std::cout << "New subgroup added" << std::endl;
      // } else {
      //   std::cout << "Subgroup already exists" << std::endl;
      // }

      if (emplace_res.second && next_node.indices.size() != group_size) {
        // Keep new node and continue the search deeper
        // std::cout << "Continuing search deeper" << std::endl;
        ++next_it;
        iters.push_back(next_it);
      } else {
        // Pop the last node, advance member iterator, continue the search
        // std::cout << "Backtracking" << std::endl;
        nodes.pop_back();
        ++(iters.back());
      }
      // std::cout << std::endl;
    }

    // Always add the trivial subgroup
    subgroups.emplace(std::set<Index>{0}, std::set<Index>{0});

    std::cout << "# subgroups found: " << subgroups.size() << std::endl;
    std::cout << "end MakeAllSubgroupsFromGenerators" << std::endl;
  }
};

/// \brief A subset of a head group, for purposes of calculating properties
///
/// By design, this class has the semantics of an immutable object; all
/// modifying operations return a new Subset object rather than modifying
/// the existing object. Some properties are lazily evaluated and cached.
template <typename ElementType>
class Subset {
 public:
  /// \brief Constructor
  ///
  /// \param group The head group
  /// \param indices Indices of elements into `group` that form the subset
  Subset(std::shared_ptr<Group<ElementType> const> const &group,
         std::set<Index> indices)
      : m_group(group),
        m_indices(std::move(indices)),
        m_is_group(std::nullopt),
        m_is_normal(std::nullopt),
        m_maximal_cyclic_generators(std::nullopt),
        m_maximal_cyclic_subgroups(std::nullopt) {
    // Validate group is not null
    if (m_group == nullptr) {
      throw std::runtime_error(
          "Error in CASM::group::Subset constructor: group is null.");
    }

    // Validate group is the head group
    if (m_group->head_group != nullptr) {
      throw std::runtime_error(
          "Error in CASM::group::Subset constructor: "
          "`group` is not the head group.");
    }

    // Validate indices are in valid range
    if (m_indices.size() > 0) {
      Index min_index = *m_indices.begin();
      if (min_index < 0) {
        throw std::runtime_error(
            "Error in CASM::group::Subset constructor: group index out of "
            "range.");
      }

      Index max_index = *m_indices.rbegin();
      if (max_index >= m_group->element.size()) {
        throw std::runtime_error(
            "Error in CASM::group::Subset constructor: group index out of "
            "range.");
      }
    }
  }

  static Subset from_generators(
      std::shared_ptr<Group<ElementType> const> const &group,
      std::set<Index> generators) {
    MakeSubgroupFromGenerators<ElementType> x(*group, generators);
    return Subset<ElementType>(group, x.indices);
  }

  /// \brief Access the indices comprising this subset
  std::set<Index> const &indices() const { return m_indices; }

  /// \brief Access the head group this subset belongs to
  std::shared_ptr<Group<ElementType> const> const &group() const {
    return m_group;
  }

  /// \brief Return true if the subset forms a group
  bool is_group() const {
    if (!m_is_group.has_value()) {
      m_is_group = this->_is_group();
    }
    return *m_is_group;
  }

  /// \brief Return true if the subset is normal
  ///
  /// A subset is normal if it is invariant under conjugation by elements of
  /// the head group.
  ///
  /// This returns true if and only if g*n*g^-1 is in the subset for all g
  /// in G and n in N, where G is the head group and N is this subset.
  bool is_normal() const {
    if (!m_is_normal.has_value()) {
      m_is_normal = this->_is_normal();
    }
    return *m_is_normal;
  }

  /// \brief Extend this subset to be the closure under group multiplication
  /// and return the result
  Subset<ElementType> close() const {
    if (!this->is_group()) {
      return this->_naive_close();
    }
    return Subset(*this);
  }

  /// \brief Extend this subset by adding indices from other
  Subset<ElementType> extend(Subset const &other) const {
    std::set<Index> new_indices(m_indices.begin(), m_indices.end());
    new_indices.insert(other.m_indices.begin(), other.m_indices.end());
    return Subset(m_group, new_indices);
  }

  /// \brief Extend and close
  Subset<ElementType> extend_and_close(Subset const &other) const {
    return extend(other).close();
  }

  /// \brief Check if this subset is a proper subset of `other`
  bool is_proper_subset_of(Subset const &other) const {
    if (m_indices.size() >= other.m_indices.size()) {
      return false;
    }

    // Both sets are ordered; use std::includes to check if `other` contains
    // all elements of `this`.
    return std::includes(other.m_indices.begin(), other.m_indices.end(),
                         m_indices.begin(), m_indices.end());
  }

  /// \brief Check if this subset is a proper subset of other or equal
  bool is_subset_of(Subset const &other) const {
    if (m_indices.size() > other.m_indices.size()) {
      return false;
    }

    // Both sets are ordered; use std::includes to check if `this` contains
    // all elements of `other`.
    return std::includes(other.m_indices.begin(), other.m_indices.end(),
                         m_indices.begin(), m_indices.end());
  }

  /// \brief Check if this subgroup is equal to other
  bool operator==(Subset const &other) const {
    return m_indices == other.m_indices;
  }

  /// \brief Check if this subgroup is not equal to other
  bool operator!=(Subset const &other) const {
    return m_indices != other.m_indices;
  }

  /// \brief Construct all cyclic subgroups generated by elements of
  /// this subset
  std::vector<Subset<ElementType>> const &all_cyclic_subgroups() const {
    if (!m_all_cyclic_subgroups.has_value()) {
      this->_make_all_cyclic_subgroups();
    }
    return *m_all_cyclic_subgroups;
  }

  /// \brief Construct the maximal cyclic subgroups generated by elements of
  /// this subset
  ///
  /// These are the cyclic subgroups generated by elements of the subset
  /// that are not a subgroup of any other cyclic subgroup of the subset.
  std::vector<Subset<ElementType>> const &maximal_cyclic_subgroups() const {
    if (!m_maximal_cyclic_subgroups.has_value()) {
      MakeMaximalCyclicSubgroups<ElementType> x(*m_group, m_indices);

      // Store results
      m_maximal_cyclic_subgroups = std::vector<Subset<ElementType>>();
      m_maximal_cyclic_generators = std::vector<Index>();

      for (auto const &cycle : x.maximal_cyclic_subgroups) {
        m_maximal_cyclic_subgroups->emplace_back(m_group, cycle);
      }
      m_maximal_cyclic_generators = x.generators;
    }
    return *m_maximal_cyclic_subgroups;
  }

  /// \brief Return generators for this subset
  ///
  /// All subset elements are included in the cyclic subgroups generated by
  /// these elements. Specifically, these are the elements that generated
  /// the maximal cyclic subgroups.
  std::vector<Index> const &maximal_cyclic_generators() const {
    if (!m_maximal_cyclic_generators.has_value()) {
      this->maximal_cyclic_subgroups();
    }
    return *m_maximal_cyclic_generators;
  }

  /// \brief Return a minimal set of generators for this subset
  ///
  /// Method:
  /// - Find the maximal cyclic subgroups of this subset.
  /// - Iteratively add the generators of those cyclic subgroups that add
  ///   the most new conjugacy classes to the generated set until all
  ///   subset elements are generated.
  ///
  /// Notes:
  /// - This is a small set of generators that generates all elements of the
  ///   subset when closed under multiplication.
  /// - It may not be the minimum set of generators.
  /// - I have not proven if this is always a minimal set of generators (i.e.,
  ///   removing any generator makes it impossible to generate all elements of
  ///   the subset), but it does generate small sets in practice.
  std::set<Index> const &minimal_generators() const {
    if (!m_minimal_generators.has_value()) {
      MakeMinimalSubsetGenerators x(*m_group, m_indices);
      m_minimal_generators = x.generators;
    }
    return *m_minimal_generators;
  }

  /// \brief Return all subgroups
  ///
  /// Uses a depth-first search for combinations of subset elements to use as
  /// subgroup generators
  ///
  /// Notes:
  /// - This subset should be a group
  ///
  std::vector<Subset<ElementType>> const &all_subgroups() const {
    if (!m_all_subgroups.has_value()) {
      MakeAllSubgroupsFromGenerators<ElementType> x(*m_group, m_indices);
      // std::cout << "# of subgroups found: " << x.subgroups.size() <<
      // std::endl; Store results
      m_all_subgroups_generators = std::vector<std::set<Index>>();
      m_all_subgroups = std::vector<Subset<ElementType>>();

      Index i = 0;
      for (auto const &res : x.subgroups) {
        m_all_subgroups_generators->push_back(res.second);
        m_all_subgroups->emplace_back(m_group, res.first);

        // std::cout << "- " << i << ": ";
        // for (Index j : res.first) {
        //   std::cout << j << " ";
        // }
        // std::cout << std::endl;
        ++i;
      }
    }
    return *m_all_subgroups;
  }

  /// \brief Return generators for all subgroups
  ///
  /// Uses a depth-first search for combinations of subset elements to use as
  /// subgroup generators
  ///
  /// Notes:
  /// - This subset should be a group
  ///
  std::vector<std::set<Index>> const &all_subgroups_generators() const {
    if (!m_all_subgroups.has_value()) {
      this->all_subgroups();
    }
    return *m_all_subgroups_generators;
  }

 private:
  bool _is_group() const {
    // Check identity
    if (m_indices.count(0) == 0) {
      return false;
    }

    // Check closure
    for (auto i : m_indices) {
      for (auto j : m_indices) {
        Index product_index = m_group->mult(i, j);
        if (m_indices.count(product_index) == 0) {
          return false;
        }
      }
    }

    // Check inverses
    for (auto i : m_indices) {
      Index inv_index = m_group->inv(i);
      if (m_indices.count(inv_index) == 0) {
        return false;
      }
    }

    return true;
  }

  /// \brief Check if this is a normal subgroup of the group
  bool _is_normal() const {
    // For each g in G and each n in N, check if g*n*g^-1 is in N
    for (auto g_index : m_group->head_group_index) {
      Index g_inv_index = m_group->inv(g_index);
      for (auto n_index : m_indices) {
        Index conjugate_index =
            m_group->mult(g_index, m_group->mult(n_index, g_inv_index));
        if (m_indices.count(conjugate_index) == 0) {
          return false;
        }
      }
    }
    return true;
  }

  /// \brief Modify this subset to be the closure under group
  /// multiplication, checking every combination of elements until no new
  /// elements are added
  Subset<ElementType> _naive_close() const {
    std::set<Index> new_indices(m_indices.begin(), m_indices.end());
    std::vector<Index> velements(m_indices.begin(), m_indices.end());
    for (Index i = 0; i < velements.size(); i++) {
      for (Index j = 0; j < velements.size(); j++) {
        Index product_index = m_group->mult(velements[i], velements[j]);
        if (new_indices.insert(product_index).second) {
          velements.push_back(product_index);
        }
      }
    }
    return Subset(m_group, new_indices);
  }

  /// \brief Use the generators of this and `other` to extend this subset
  /// and close it under multiplication
  ///
  /// This approach does not use `_naive_close`. The method is:
  /// - Create the union of the elements in this and other.
  /// - Create the union of the the minimal generators for this and other.
  /// - Use a queue to iteratively multiply generators and the newly added
  ///   elements until no new elements are added.
  void _extend_and_close(Subset const &other);

  void _make_all_cyclic_subgroups() const {
    m_all_cyclic_subgroups = std::vector<Subset<ElementType>>();

    for (Index i : m_indices) {
      std::set<Index> indices;
      Index prod = i;
      while (indices.insert(prod).second) {
        prod = m_group->mult(i, prod);
      }
      m_all_cyclic_subgroups->emplace_back(m_group, indices);
    }
  }

  std::shared_ptr<Group<ElementType> const> m_group;

  std::set<Index> m_indices;

  /// \brief Stores whether the subset is a group, if known
  mutable std::optional<bool> m_is_group;

  /// \brief Stores whether the subset is normal (invariant to conjugation),
  /// if known
  mutable std::optional<bool> m_is_normal;

  /// \brief A vector of all cyclic subgroups generated by elements of
  /// this subset, if known
  mutable std::optional<std::vector<Subset<ElementType>>>
      m_all_cyclic_subgroups;

  /// \brief A vector of generators for this subgroup, if known.
  ///
  /// - All subgroup elements are included in the cyclic subgroups generated
  /// by these generators
  /// - m_maximal_cyclic_generators[i] generates m_maximal_cyclic_subgroups[i]
  mutable std::optional<std::vector<Index>> m_maximal_cyclic_generators;

  /// \brief The cyclic subgroups generated by the generators
  ///
  /// - Does not include cyclic subgroups that are a subgroup of another
  ///   cyclic subgroup (i.e. the cyclic subgroup of a 4-fold rotation is
  ///   included, but the cyclic subgroup of the associated 2-fold rotation
  ///   is not included separately).
  mutable std::optional<std::vector<Subset<ElementType>>>
      m_maximal_cyclic_subgroups;

  /// \brief A minimal set of generators
  ///
  /// - This may not be the minimum set of generators
  mutable std::optional<std::set<Index>> m_minimal_generators;

  /// \brief Generators for `m_all_subgroups`
  mutable std::optional<std::vector<std::set<Index>>>
      m_all_subgroups_generators;

  /// \brief A vector of all subgroups of this subset
  mutable std::optional<std::vector<Subset<ElementType>>> m_all_subgroups;
};

}  // namespace group
}  // namespace CASM
// --- Implementation ---

#include <numeric>

namespace CASM {
namespace group {

namespace subgroups_impl {

typedef std::set<Index> CosetIndices;

/// \brief Return the unique left cosets of a subgroup
///
/// - If subgroup B of group G contains elements: (E, B1, B2, …, Bg),
/// the "left coset” of X is (X*E, X*B1, X*B2, …, X*Bg),
/// where X is an element of G.
/// - A coset need not be a subgroup.
/// - If X is an element of B, then the coset will be a subgroup of B.
/// - Two left cosets of a given subgroup either contain exactly the same
/// elements, or have no elements in common.
template <typename ElementType>
std::set<CosetIndices> _make_left_cosets(
    Group<ElementType> const &group, SubgroupIndices const &subgroup_indices) {
  std::set<CosetIndices> left_cosets;

  // each group element is only included in one coset
  std::vector<bool> check(group.element.size(), false);
  Index product_index;
  Index group_element_index = 0;
  while (group_element_index < group.element.size()) {
    if (check[group_element_index]) {
      ++group_element_index;
      continue;
    }
    CosetIndices left_coset;
    for (auto subgroup_element_index : subgroup_indices) {
      product_index = group.mult(group_element_index, subgroup_element_index);
      left_coset.insert(product_index);
      check[product_index] = true;
    }
    left_cosets.insert(left_coset);
    ++group_element_index;
  }
  return left_cosets;
}

template <typename ElementType>
SubgroupOrbit _make_subgroup_orbit(Group<ElementType> const &group,
                                   SubgroupIndices const &subgroup) {
  SubgroupOrbit orbit;
  std::set<subgroups_impl::CosetIndices> left_cosets =
      subgroups_impl::_make_left_cosets(group, subgroup);
  for (auto const &coset : left_cosets) {
    Index X_index = *coset.begin();
    Index X_inv_index = group.inv(X_index);
    SubgroupIndices equiv_subgroup;
    for (auto const &A_index : subgroup) {
      equiv_subgroup.insert(
          group.mult(X_index, group.mult(A_index, X_inv_index)));
    }
    orbit.insert(equiv_subgroup);
  }
  return orbit;
}

inline std::function<bool(SubgroupIndices const &)> _make_subgroup_count(
    std::set<SubgroupOrbit> const &subgroups) {
  return [&](SubgroupIndices const &subgroup) {
    for (auto const &orbit : subgroups) {
      if (orbit.count(subgroup)) {
        return true;
      }
    }
    return false;
  };
}

template <typename ElementType>
std::function<void(SubgroupIndices &)> _make_close_subgroup(
    Group<ElementType> const &group) {
  return [&](SubgroupIndices &subgroup) {
    std::vector<Index> vgroup(subgroup.begin(), subgroup.end());
    for (Index i = 0; i < vgroup.size(); i++) {
      for (Index j = 0; j < vgroup.size(); j++) {
        Index product_index = group.mult(vgroup[i], vgroup[j]);
        if (subgroup.insert(product_index).second) {
          vgroup.push_back(product_index);
        }
      }
    }
  };
}

}  // namespace subgroups_impl

/// \brief Return all cyclic subgroups
///
/// - A cyclic subgroup, A, is the subgroup generated by repeated
/// multiplication of a single element, i.e. A={a, a^2, a^3, ..., a^k}, where
/// a^k=E, the identity element
/// - An equivalent subgroup is {X*a*X^-1, X*(a^2)*X^-1, X*(a^3)*X^-1, ...,
/// E}, where X is an element in a coset of A
/// - The orbit of a cyclic subgroup is all distinct equivalent subgroups
///
/// \param group The group to find cyclic subgroups of
/// \returns A set of orbits of cyclic subgroups
///
template <typename ElementType>
std::set<SubgroupOrbit> make_cyclic_subgroups(Group<ElementType> const &group) {
  using namespace subgroups_impl;

  std::set<SubgroupOrbit> cyclic_subgroups;
  Index group_element_index = 0;
  while (group_element_index < group.element.size()) {
    // Make cyclic subgroup of element `group_element_index`
    SubgroupIndices cyclic_subgroup;
    cyclic_subgroup.insert(group_element_index);
    Index product_index = group_element_index;
    while (product_index != 0) {
      product_index = group.mult(group_element_index, product_index);
      cyclic_subgroup.insert(product_index);
    }

    // Make orbit of subgroups equivalent to `cyclic_subgroup` && Insert orbit
    cyclic_subgroups.insert(_make_subgroup_orbit(group, cyclic_subgroup));

    ++group_element_index;
  }
  return cyclic_subgroups;
}

// template <typename ElementType>
// struct LatticeElement {
//   // Constructor
//   LatticeElement(SubgroupIndices const &_subgroup) : subgroup(_subgroup) {}
//
//   std::shared_ptr<Group<ElementType> const> head_group;
//   SubgroupIndices subgroup;
//   std::vector<std::shared_ptr<LatticeElement>> subgroups;
//   std::vector<std::shared_ptr<LatticeElement>> parent_groups;
//
//   LatticeElement extend(LatticeElement const &other) {
//     SubgroupIndices new_subgroup = this->subgroup;
//     new_subgroup.insert(other.subgroup.begin(), other.subgroup.end());
//     return LatticeElement(new_subgroup);
//   }
// };
//
// template <typename ElementType>
// struct Lattice {
//   // Constructor
//   Lattice(std::shared_ptr<Group<ElementType> const> const &_head_group)
//       : head_group(_head_group) {}
//
//   // Constructor
//   Lattice(std::shared_ptr<Group<ElementType> const> const &_head_group,
//           std::vector<std::shared_ptr<LatticeElement<ElementType>>> const
//               &_elements)
//       : head_group(_head_group), elements(_elements) {}
//
//   std::shared_ptr<Group<ElementType> const> head_group;
//
//   std::vector<std::shared_ptr<LatticeElement<ElementType>>> elements;
// };
//
// template <typename ElementType>
// Lattice make_subgroup_lattice(Group<ElementType> const &group,
//                               SubgroupIndices const &subgroup) {
//   // Make cyclic subgroups of the subgroup
//   using namespace subgroups_impl;
//
//   // Validate that `group` element 0 is the identity element
//   if (group.inv(0) != 0 || group.mult(0, 0) != 0) {
//     throw std::runtime_error(
//         "Error in make_subgroup_lattice: group element 0 is not the
//         identity " "element.");
//   }
//
//   SubgroupIndices identity_subgroup = {0};
//   std::set<SubgroupIndices> cyclic_subgroups;
//   for (Index generator : subgroup) {
//     // Make cyclic subgroup of element `generator`
//     SubgroupIndices cyclic_subgroup;
//     cyclic_subgroup.insert(generator);
//     Index product_index = generator;
//     while (product_index != 0) {
//       product_index = group.mult(generator, product_index);
//       cyclic_subgroup.insert(product_index);
//     }
//
//     // If cyclic_subgroup is the trivial subgroup, set identity_subgroup
//     if (cyclic_subgroup.size() == 1 && *(cyclic_subgroup.begin()) == 0) {
//       identity_subgroup = cyclic_subgroup;
//       continue;
//     }
//
//     // Insert subgroup
//     cyclic_subgroups.insert(cyclic_subgroup);
//   }
//
//   Lattice lattice;
//   lattice.elements.push_back(
//       std::make_shared<LatticeElement>(identity_subgroup));
// }

/// \brief Return all subgroups
///
/// Method:
/// - Start with m_subgroups = m_small_subgroups, then add new subgroups by
/// finding the closure of a union of a large_group and a small_group.
/// - If the the new large_group is unique, add it as a large_group.
/// - Repeat for all (large_group, small_group) pairs, until no new
/// m_subgroups are found.
///
/// Note:
/// - This is probably not the fastest algorithm, but it is complete
///
/// \param group The group to find subgroups of
/// \returns A set of orbits of all subgroups
///
template <typename ElementType>
std::set<SubgroupOrbit> make_all_subgroups(Group<ElementType> const &group) {
  using namespace subgroups_impl;
  std::set<SubgroupOrbit> small_subgroups = make_cyclic_subgroups(group);
  std::set<SubgroupOrbit> all_subgroups = small_subgroups;

  // functor to close incomplete subgroup
  auto _close_subgroup = _make_close_subgroup(group);

  // functor to find if any orbit contains a particular subgroup
  auto _subgroup_count = _make_subgroup_count(all_subgroups);

  auto all_subgroups_it = all_subgroups.begin();
  while (all_subgroups_it != all_subgroups.end()) {
    for (auto const &small_subgroups_orbit : small_subgroups) {
      for (auto const &small_subgroups_equiv : small_subgroups_orbit) {
        // Combine an existing subgroup and a small (cyclic) subgroup
        SubgroupIndices subgroup = *(all_subgroups_it->begin());
        Index init_size = subgroup.size();
        subgroup.insert(small_subgroups_equiv.begin(),
                        small_subgroups_equiv.end());
        if (subgroup.size() == init_size) continue;

        // Find group closure
        _close_subgroup(subgroup);

        // If subgroup already exists in all_subgroups, continue
        if (_subgroup_count(subgroup)) continue;

        // Else, make orbit and insert
        all_subgroups.insert(_make_subgroup_orbit(group, subgroup));
      }
    }
    ++all_subgroups_it;
  }
  return all_subgroups;
}

/// A functor which makes all cyclic subgroups at construction and then
/// returns them when called:
template <typename ElementType>
class MakeCyclicSubgroups {
 public:
  MakeCyclicSubgroups(std::shared_ptr<Group<ElementType> const> group)
      : m_group(group),
        m_subgroups_constructed(false),
        m_cyclic_subgroups(std::make_shared<std::set<SubgroupOrbit>>()) {}

  std::set<SubgroupOrbit> operator()() const {
    if (!m_subgroups_constructed) {
      *m_cyclic_subgroups = make_cyclic_subgroups(*m_group);
      m_subgroups_constructed = true;
    }
    return *m_cyclic_subgroups;
  }

 private:
  std::shared_ptr<Group<ElementType> const> m_group;

  mutable bool m_subgroups_constructed;

  mutable std::shared_ptr<std::set<SubgroupOrbit>> m_cyclic_subgroups;
};

//// A functor which makes all subgroups at construction and then
/// returns them when called:
template <typename ElementType>
class MakeAllSubgroups {
 public:
  MakeAllSubgroups(std::shared_ptr<Group<ElementType> const> group)
      : m_group(group),
        m_subgroups_constructed(std::make_shared<bool>(false)),
        m_all_subgroups(std::make_shared<std::set<SubgroupOrbit>>()) {}

  std::set<SubgroupOrbit> operator()() const {
    if (*m_subgroups_constructed == false) {
      *m_all_subgroups = make_all_subgroups(*m_group);
      *m_subgroups_constructed = true;
    }
    return *m_all_subgroups;
  }

 private:
  std::shared_ptr<Group<ElementType> const> m_group;

  mutable std::shared_ptr<bool> m_subgroups_constructed;

  mutable std::shared_ptr<std::set<SubgroupOrbit>> m_all_subgroups;
};

/// \brief Make the invariant subgroup for each orbit element, as
///     indices of group elements
///
/// \param equivalence_map The indices equivalence_map[i] are the
///     indices of the group elements which map orbit element 0 onto
///     orbit element i
/// \param group The group used to generate the orbit
///
/// \returns invariant_subgroups, The indices invariant_subgroups[i] are
///     the indices of the group elements which leave orbit element i
///     invariant.
template <typename GroupElementType>
std::vector<SubgroupIndices> make_invariant_subgroups(
    std::vector<std::vector<Index>> const &equivalence_map,
    Group<GroupElementType> const &group) {
  /// The first row of equivalence_map is the invariant subgroup
  /// of the first element in the orbit
  ///
  /// Invariant subgroups of subsequent orbit elements are constructed
  /// using:
  ///   inv_subgrp(i) = eq_map(i,0) * eq_map(0,j) * inverse(eq_map(i,0)), for
  ///   all j
  ///
  /// inverse(eq_map(i,0)): transforms orbit element i back to orbit element 0
  /// eq_map(0,j): invariant transformation of orbit element 0
  /// eq_map(i,0): transforms orbit element 0 to orbit element i
  ///
  std::vector<SubgroupIndices> invariant_subgroups;
  if (!equivalence_map.size()) {
    return invariant_subgroups;
  }

  // first row of equivalence map is an invariant subgroup
  {
    SubgroupIndices subgroup;
    for (Index e_0j : equivalence_map[0]) {
      subgroup.insert(e_0j);
    }
    invariant_subgroups.emplace_back(std::move(subgroup));
  }

  // first column can be used to generate others
  for (Index i = 1; i < equivalence_map.size(); ++i) {
    if (!equivalence_map[i].size()) {
      throw std::runtime_error(
          "Error in make_invariant_subgroups: failed due to empty row in "
          "equivalence_map");
    }
    SubgroupIndices subgroup;
    Index e_i0 = equivalence_map[i][0];
    for (Index e_0j : equivalence_map[0]) {
      subgroup.insert(group.mult(e_i0, group.mult(e_0j, group.inv(e_i0))));
    }
    invariant_subgroups.emplace_back(std::move(subgroup));
  }
  return invariant_subgroups;
}

}  // namespace group
}  // namespace CASM

#endif
