#ifndef CASM_group_Group
#define CASM_group_Group

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <set>

#include "casm/configuration/group/definitions.hh"
#include "casm/global/threads.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace group {

/// \brief Holds multiplication table and derived data
///
/// Notes:
/// - The constructor generates the inverse_index and class_index data from the
///   multiplication table.
/// - When the inverse indices are generated, a check is performed to ensure
///   that the multiplication table is valid. This validates that:
///   - the table is square
///   - the first element is the identity element
///   - each element has an inverse
///   - the table is closed under multiplication
///   - each row and column contains each element exactly once
///
struct GenericGroup {
  /// \brief Construct a head group
  GenericGroup(MultiplicationTable _multiplication_table);

  /// \brief Construct a subgroup
  GenericGroup(std::shared_ptr<GenericGroup const> _head_group_ptr,
               std::set<Index> _head_group_index);

  /// \brief Get the pointer to the head group
  std::shared_ptr<GenericGroup const> head() const { return head_group_ptr; }

  /// \brief If this is a subgroup, indicates the head group; if this is a head
  /// group, then this is empty
  ///
  /// Notes:
  /// - This is a copy of the pointer Group<ElementType>::head_group when this
  ///   is used as the base class of Group<ElementType>. The duplication is to
  ///   avoid breaking existing code.
  std::shared_ptr<GenericGroup const> const head_group_ptr;

  /// \brief Specifies the head group index for each element (guaranteed sorted)
  ///
  /// If this is the head group, then:
  ///
  ///     this->head_group_index = [0, 1, 2, ...]
  ///
  /// If this is a sub group, then in a derived Group<ElementType>:
  ///
  ///     this->element[i] == head()->element[this->head_group_index[i]]
  ///
  /// Or, for example, for a subgroup of a factor group:
  ///
  ///     this->element[i] == <translation> *
  ///     head()->element[this->head_group_index[i]]
  ///
  std::vector<Index> const head_group_index;

  /// \brief Specifies the multiplication table for the elements
  ///
  /// element[k] == element[i] * element[j],
  /// where k = multiplication_table[i][j]
  MultiplicationTable const multiplication_table;

  /// \brief Specifies the index of the inverse element
  ///
  ///     I == element[i] * element[inverse_index[i]]
  ///       == element[inverse_index[i]] * element[i]
  std::vector<Index> const inverse_index;

  /// \brief The conjugacy class of each element
  ///
  /// The `i`-th element is in the `cc`-th class, where `cc = class_index[i]`.
  std::vector<Index> class_index;

  std::size_t size() const { return head_group_index.size(); }

  /// \brief Use the multiplication table
  ///
  /// \param i,j Element indices
  /// \returns k, where element[k] == element[i] * element[j]
  Index mult(Index i, Index j) const { return multiplication_table[i][j]; }

  /// \brief Get the inverse element index
  ///
  /// \param i Element index
  /// \returns i_inv, The index of the inverse element of element i
  Index inv(Index i) const { return inverse_index[i]; }

  /// \brief Get the conjugacy class index of an element
  ///
  /// \param i Element index
  /// \return cc, The index of the conjugacy class containing element i
  Index class_of(Index i) const { return class_index[i]; }

  bool is_subgroup() const { return head_group_ptr != nullptr; }
};

/// \brief Holds group elements and multiplication table
template <typename ElementType>
struct Group : public GenericGroup {
  /// \brief Construct a head group
  Group(std::vector<ElementType> const &_element,
        MultiplicationTable const &_multiplication_table);

  /// \brief Construct a subgroup
  Group(std::shared_ptr<Group const> const &_head_group,
        std::set<Index> const &_head_group_index);

  /// \brief Construct a subgroup
  Group(std::shared_ptr<Group const> const &_head_group,
        std::vector<ElementType> const &_element,
        std::set<Index> const &_head_group_index);

  /// \brief Get the pointer to the head group
  ///
  /// This version casts from GenericGroup::head_group_ptr to the
  /// appropriate type.
  std::shared_ptr<Group const> head() const {
    return std::static_pointer_cast<Group const>(head_group_ptr);
  }

  /// \brief If this is a subgroup, indicates the head group; if this is a head
  /// group, then this is empty.
  ///
  /// \deprecated Use head() instead.
  std::shared_ptr<Group const> const head_group;

  /// \brief Specifies the group elements
  std::vector<ElementType> const element;
};

template <typename ElementType,
          typename MultiplyFunctionType = std::multiplies<ElementType>,
          typename EqualToFunctionType = std::equal_to<ElementType>>
Group<ElementType> make_group(
    std::vector<ElementType> const &element,
    MultiplyFunctionType multiply_f = MultiplyFunctionType(),
    EqualToFunctionType equal_to_f = EqualToFunctionType(),
    bool sort_by_class = false);

/// \brief Determine conjugacy classes
std::vector<std::vector<Index>> make_conjugacy_classes(
    GenericGroup const &group);

/// \brief Make map of element index to conjugacy class index
std::vector<Index> make_element_to_class(GenericGroup const &group);

}  // namespace group
}  // namespace CASM

// --- Implementation ---

#include <numeric>

namespace CASM {
namespace group {

namespace Group_impl {

inline std::vector<Index> _identity_indices(Index n) {
  std::vector<Index> result(n);
  std::iota(result.begin(), result.end(), 0);
  return result;
}

template <typename ElementType>
std::vector<ElementType> _make_subgroup_elements(
    std::shared_ptr<Group<ElementType> const> const &_head_group,
    std::set<Index> const &_head_group_index) {
  std::vector<ElementType> result;
  for (Index index : _head_group_index) {
    result.push_back(_head_group->element[index]);
  }
  return result;
}

inline MultiplicationTable _make_subgroup_multiplication_table(
    std::shared_ptr<GenericGroup const> const &_head_group,
    std::set<Index> const &_head_group_index) {
  MultiplicationTable result(_head_group_index.size());
  MultiplicationTable const &head_group_table =
      _head_group->multiplication_table;
  Index N = head_group_table.size();

  for (Index index : _head_group_index) {
    if (index >= N) {
      throw std::runtime_error(
          "Error in Group constructor: head group index >= head group "
          "multiplication table size");
    }
  }

  Index row = 0;
  for (Index i = 0; i < N; ++i) {
    if (head_group_table[i].size() != N) {
      throw std::runtime_error(
          "Error in Group constructor: head group multiplication table is not "
          "square");
    }
    if (!_head_group_index.count(i)) {
      continue;
    }

    for (Index j = 0; j < N; ++j) {
      if (!_head_group_index.count(j)) {
        continue;
      }
      auto it = _head_group_index.find(head_group_table[i][j]);
      if (it == _head_group_index.end()) {
        throw std::runtime_error(
            "Error in Group constructor: subgroup is not closed according to "
            "the head group multiplication table.");
      }
      Index subgroup_entry = std::distance(_head_group_index.begin(), it);
      result[row].push_back(subgroup_entry);
    }

    ++row;
  }
  return result;
}

/// \brief Collect indices of the inverse elements in a group using the
/// multiplication table
///
/// Notes:
/// - This also validates that the multiplication table is square, the
///   first element is the identity, that each element has an inverse, and that
///   each element appears exactly once in each row and column.
/// - Requires that identity element corresponds to index 0
inline std::vector<Index> _make_inverse_index(
    MultiplicationTable const &multiplication_table) {
  std::vector<Index> index_inverse;

  Index N = multiplication_table.size();
  for (Index i = 0; i < N; ++i) {
    if (multiplication_table[i].size() != N) {
      throw std::runtime_error(
          "Error in Group constructor: multiplication table is not square");
    }
  }

  // Check that element 0 is identity. This requires that
  // multiplication_table[0][i] == i and multiplication_table[i][0] == i
  for (Index i = 0; i < N; ++i) {
    if (multiplication_table[0][i] != i) {
      std::cout << "Multiplicaiton table row 0: ";
      for (Index j = 0; j < N; ++j) {
        std::cout << multiplication_table[0][j] << " ";
      }
      std::cout << std::endl;
      throw std::runtime_error(
          "Error in make_inverse_index: multiplication table identity error");
    }
    if (multiplication_table[i][0] != i) {
      std::cout << "Multiplicaiton table column 0: ";
      for (Index j = 0; j < N; ++j) {
        std::cout << multiplication_table[j][0] << " ";
      }
      std::cout << std::endl;
      throw std::runtime_error(
          "Error in make_inverse_index: multiplication table identity error");
    }
  }

  // Find inverse elements
  for (auto const &row : multiplication_table) {
    auto begin = std::begin(row);
    auto end = std::end(row);
    auto it = std::find(begin, end, 0);
    if (it == end) {
      throw std::runtime_error(
          "Error in make_inverse_index: no inverse element");
    }
    index_inverse.push_back(std::distance(begin, it));
  }

  // Validate that each element appears exactly once in each row and column
  std::vector<bool> row_check(N, false);
  std::vector<bool> col_check(N, false);
  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      Index row_entry = multiplication_table[i][j];
      Index col_entry = multiplication_table[j][i];
      if (row_entry < 0 || row_entry >= N) {
        throw std::runtime_error(
            "Error in make_inverse_index: multiplication table entry out of "
            "range");
      }
      if (col_entry < 0 || col_entry >= N) {
        throw std::runtime_error(
            "Error in make_inverse_index: multiplication table entry out of "
            "range");
      }
      if (row_check[row_entry]) {
        throw std::runtime_error(
            "Error in make_inverse_index: duplicate entry in multiplication "
            "table row");
      }
      if (col_check[col_entry]) {
        throw std::runtime_error(
            "Error in make_inverse_index: duplicate entry in multiplication "
            "table column");
      }
      row_check[row_entry] = true;
      col_check[col_entry] = true;
    }
    std::fill(row_check.begin(), row_check.end(), false);
    std::fill(col_check.begin(), col_check.end(), false);
  }

  return index_inverse;
}

}  // namespace Group_impl

/// \brief Construct a head group
inline GenericGroup::GenericGroup(MultiplicationTable _multiplication_table)
    : head_group_ptr(nullptr),
      head_group_index(
          Group_impl::_identity_indices(_multiplication_table.size())),
      multiplication_table(_multiplication_table),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {
  class_index = make_element_to_class(*this);
}

/// \brief Construct a subgroup
inline GenericGroup::GenericGroup(
    std::shared_ptr<GenericGroup const> _head_group_ptr,
    std::set<Index> _head_group_index)
    : head_group_ptr(std::move(_head_group_ptr)),
      head_group_index(_head_group_index.begin(), _head_group_index.end()),
      multiplication_table(Group_impl::_make_subgroup_multiplication_table(
          head_group_ptr, _head_group_index)),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {
  class_index = make_element_to_class(*this);
}

/// \brief Construct a head group
///
/// \params _element Group elements, expected to be closed and sorted as desired
/// \params _multiplication_table Contains indices of products,
///     `_multiplication_table[i][j] == _element[i] * _element[j]`.
///
/// Notes:
/// - Use `make_group` to build the multiplication table from known
///   multiplication and equals_to operations.
template <typename ElementType>
Group<ElementType>::Group(std::vector<ElementType> const &_element,
                          MultiplicationTable const &_multiplication_table)
    : GenericGroup(_multiplication_table),
      head_group(nullptr),
      element(_element) {}

/// \brief Construct a subgroup
///
/// \params _head_group The group that is the head group of this subgroup.
/// \params _head_group_index Contains indices into `_head_group->element` of
///     the members of the subgroup.
///
template <typename ElementType>
Group<ElementType>::Group(
    std::shared_ptr<Group<ElementType> const> const &_head_group,
    std::set<Index> const &_head_group_index)
    : GenericGroup(_head_group, _head_group_index),
      head_group(_head_group),
      element(Group_impl::_make_subgroup_elements(_head_group,
                                                  _head_group_index)) {}

/// \brief Construct a subgroup
///
/// Note:
/// This constructor allows representing subgroups of the space group with
/// reference to the factor group. For example, cluster invariant
/// groups need SymOp that have the correct translation to leave the
/// cluster invariant, which may have a different translation than
/// the corresponding element in the factor group elements list.
///
/// \params _head_group The group that is the head group of this subgroup.
///     the members of the subgroup.
/// \params _element Group elements, expected to be closed and in order
///     consistent with _head_group_index.
/// \params _head_group_index Contains indices into `_head_group->element` of
///     the members of the subgroup.
///
template <typename ElementType>
Group<ElementType>::Group(
    std::shared_ptr<Group<ElementType> const> const &_head_group,
    std::vector<ElementType> const &_element,
    std::set<Index> const &_head_group_index)
    : GenericGroup(_head_group, _head_group_index),
      head_group(_head_group),
      element(_element) {}

/// \brief Construct a head group with optional sorting by class
///
/// Notes:
/// - This builds the multiplication table from the provided multiplication and
///   equality functions.
/// - If `sort_by_class` is true, then the elements are sorted by conjugacy
///   class, using the initial order to break ties within each class and to
///   order the classes.
///
/// \param element The group elements, expected to be closed and with the first
///     element being identity.
/// \param multiply_f A function that takes two elements and returns their
///     product.
/// \param equal_to_f A function that takes two elements and returns true if
///     they are equal, and false otherwise.
/// \param sort_by_class If true, then the elements are sorted by conjugacy
///     class, using the initial order to break ties within each class and
///     to order the classes.
///
/// \returns The Group object.
///
template <typename ElementType, typename MultiplyFunctionType,
          typename EqualToFunctionType>
Group<ElementType> make_group(std::vector<ElementType> const &element,
                              MultiplyFunctionType multiply_f,
                              EqualToFunctionType equal_to_f,
                              bool sort_by_class) {
  Index size = element.size();

  // multi-threaded version:

  // preallocate a square table so each thread can safely write to distinct rows
  MultiplicationTable multiplication_table(size, std::vector<Index>(size));

  auto worker = [&](Index start, Index end, Index thread_id) {
    for (Index i = start; i < end; ++i) {
      for (Index j = 0; j < size; ++j) {
        ElementType product = multiply_f(element[i], element[j]);
        auto it = std::find_if(
            element.begin(), element.end(),
            [&](ElementType const &lhs) { return equal_to_f(lhs, product); });
        if (it == element.end()) {
          request_stop();
          throw std::runtime_error(
              "Error in CASM::group::make_group: Failed to construct "
              "multiplication table");
        }
        multiplication_table[i][j] =
            static_cast<Index>(std::distance(element.begin(), it));
      }
    }
  };

  threaded_run(size, worker);

  if (!sort_by_class) {
    return Group<ElementType>(element, multiplication_table);
  }

  /// Get conjugacy classes to sort elements by class
  auto tmp = Group<ElementType>(element, multiplication_table);
  auto conjugacy_classes = make_conjugacy_classes(tmp);

  /// Build the sorted elements and a lookup from input index to sorted index
  std::vector<ElementType> sorted_element;
  std::vector<Index> input_index_to_sorted_index(size);
  for (auto const &cclass : conjugacy_classes) {
    for (Index index : cclass) {
      input_index_to_sorted_index[index] = sorted_element.size();
      sorted_element.push_back(element[index]);
    }
  }

  /// Build the multiplication table for the sorted elements from the
  /// multiplication table for the elements in the initial order
  MultiplicationTable sorted_table(size, std::vector<Index>(size));
  auto f = [&](Index i) { return input_index_to_sorted_index[i]; };
  for (Index i = 0; i < size; ++i) {
    for (Index j = 0; j < size; ++j) {
      sorted_table[f(i)][f(j)] = f(multiplication_table[i][j]);
    }
  }

  return Group<ElementType>(sorted_element, sorted_table);
}

/// \brief Determine conjugacy classes
///
/// Notes:
/// - This maintains the order of elements in a class according to their
///   original order in the group.
/// - Classes are ordered according to the index of the first element in each
///   class.
///
/// \returns conjugacy_classes, where conjugacy_classes[i] is a vector of
///     the indices of elements in class 'i'
///
inline std::vector<std::vector<Index>> make_conjugacy_classes(
    GenericGroup const &group) {
  std::vector<std::vector<Index>> conjugacy_classes;

  // check if operation i is in an existing class
  auto is_in_existing_class = [&](Index i) {
    for (auto const &cclass : conjugacy_classes) {
      if (contains(cclass, i)) {
        return true;
      }
    }
    return false;
  };

  Index group_size = group.multiplication_table.size();
  for (Index i = 0; i < group_size; i++) {
    if (is_in_existing_class(i)) continue;

    std::set<Index> curr_class;
    for (Index j = 0; j < group_size; j++) {
      curr_class.insert(group.mult(j, group.mult(i, group.inv(j))));
    }
    conjugacy_classes.emplace_back(curr_class.begin(), curr_class.end());
  }

  return conjugacy_classes;
}

/// \brief Make map of element index to conjugacy class index
inline std::vector<Index> make_element_to_class(GenericGroup const &group) {
  std::vector<std::vector<Index>> conjugacy_classes =
      make_conjugacy_classes(group);
  std::vector<Index> element_to_class(group.multiplication_table.size());
  for (Index class_index = 0; class_index < conjugacy_classes.size();
       ++class_index) {
    for (Index element_index : conjugacy_classes[class_index]) {
      element_to_class[element_index] = class_index;
    }
  }
  return element_to_class;
}

}  // namespace group
}  // namespace CASM

#endif
