#ifndef CASM_basic_symmetry_Group
#define CASM_basic_symmetry_Group

#include <memory>
#include <set>

#include "casm/basic_symmetry/definitions.hh"

namespace CASM {
namespace basic_symmetry {

/// \brief Holds group elements and multiplication table
template <typename ElementType>
struct Group {
  /// \brief Construct a head group
  Group(std::vector<ElementType> const &_element,
        MultiplicationTable const &_multiplication_table);

  /// \brief Construct a subgroup
  Group(std::shared_ptr<Group const> const &_head_group,
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
};

template <typename ElementType,
          typename MultiplyFunctionType = std::multiplies<ElementType>,
          typename EqualToFunctionType = std::equal_to<ElementType>>
Group<ElementType> make_group(
    std::vector<ElementType> const &element,
    MultiplyFunctionType multiply_f = MultiplyFunctionType(),
    EqualToFunctionType equal_to_f = EqualToFunctionType());

}  // namespace basic_symmetry
}  // namespace CASM

// --- Implementation ---

#include <numeric>

namespace CASM {
namespace basic_symmetry {

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
            "Error in CASM::basic_symmetry::make_group: Failed to construct "
            "multiplication table");
      }

      multiplication_table[i].push_back(std::distance(begin, it));
    }
  }
  return Group<ElementType>(element, multiplication_table);
}

}  // namespace basic_symmetry
}  // namespace CASM

#endif
