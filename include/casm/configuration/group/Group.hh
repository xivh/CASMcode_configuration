#ifndef CASM_group_Group
#define CASM_group_Group

#include <memory>
#include <set>

#include "casm/configuration/group/definitions.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {
namespace group {

/// \brief Holds group elements and multiplication table
template <typename ElementType>
struct Group {
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

  /// \brief If this is a subgroup, indicates the head group; if this is a head
  /// group, then this is empty
  std::shared_ptr<Group const> const head_group;

  /// \brief Specifies the group elements
  std::vector<ElementType> const element;

  /// \brief Specifies the head group index for each element (guaranteed sorted)
  ///
  /// If this is the head group, then:
  ///
  ///     this->head_group_index = [0, 1, 2, ...]
  ///
  /// If this is a sub group, then:
  ///
  ///     this->element[i] == head_group->element[this->head_group_index[i]]
  ///
  /// Or, for example, for a subgroup of a factor group:
  ///
  ///     this->element[i] == <translation> *
  ///     head_group->element[this->head_group_index[i]]
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
};

template <typename ElementType,
          typename MultiplyFunctionType = std::multiplies<ElementType>,
          typename EqualToFunctionType = std::equal_to<ElementType>>
Group<ElementType> make_group(
    std::vector<ElementType> const &element,
    MultiplyFunctionType multiply_f = MultiplyFunctionType(),
    EqualToFunctionType equal_to_f = EqualToFunctionType());

/// \brief Determine conjugacy classes
template <typename ElementType>
std::vector<std::vector<Index>> make_conjugacy_classes(
    Group<ElementType> const &group);

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

template <typename ElementType>
MultiplicationTable _make_subgroup_multiplication_table(
    std::shared_ptr<Group<ElementType> const> const &_head_group,
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
  return index_inverse;
}

}  // namespace Group_impl

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
    : head_group(nullptr),
      element(_element),
      head_group_index(Group_impl::_identity_indices(element.size())),
      multiplication_table(_multiplication_table),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {}

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
    : head_group(_head_group),
      element(
          Group_impl::_make_subgroup_elements(_head_group, _head_group_index)),
      head_group_index(_head_group_index.begin(), _head_group_index.end()),
      multiplication_table(Group_impl::_make_subgroup_multiplication_table(
          _head_group, _head_group_index)),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {}

/// \brief Construct a subgroup
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
    : head_group(_head_group),
      element(_element),
      head_group_index(_head_group_index.begin(), _head_group_index.end()),
      multiplication_table(Group_impl::_make_subgroup_multiplication_table(
          _head_group, _head_group_index)),
      inverse_index(Group_impl::_make_inverse_index(multiplication_table)) {}

template <typename ElementType, typename MultiplyFunctionType,
          typename EqualToFunctionType>
Group<ElementType> make_group(std::vector<ElementType> const &element,
                              MultiplyFunctionType multiply_f,
                              EqualToFunctionType equal_to_f) {
  Index size = element.size();
  MultiplicationTable multiplication_table(size);
  auto begin = element.begin();
  auto end = element.end();
  for (Index i = 0; i < size; ++i) {
    for (Index j = 0; j < size; ++j) {
      ElementType product = multiply_f(element[i], element[j]);
      auto unary_f = [&](ElementType const &lhs) {
        return equal_to_f(lhs, product);
      };
      auto it = std::find_if(begin, end, unary_f);
      if (it == end) {
        throw std::runtime_error(
            "Error in CASM::group::make_group: Failed to construct "
            "multiplication table");
      }

      multiplication_table[i].push_back(std::distance(begin, it));
    }
  }
  return Group<ElementType>(element, multiplication_table);
}

/// \brief Determine conjugacy classes
///
/// \returns conjugacy_classes, where conjugacy_classes[i] is a vector of
///     the indices of elements in class 'i'
///
template <typename ElementType>
std::vector<std::vector<Index>> make_conjugacy_classes(
    Group<ElementType> const &group) {
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

  Index group_size = group.element.size();
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

}  // namespace group
}  // namespace CASM

#endif
