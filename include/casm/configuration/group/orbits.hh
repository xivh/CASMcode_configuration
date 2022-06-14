#ifndef CASM_group_orbits
#define CASM_group_orbits

#include <memory>
#include <set>

#include "casm/configuration/group/definitions.hh"

namespace CASM {
namespace group {

template <typename OrbitElementType, typename GroupElementIt,
          typename CompareType, typename CopyApplyType>
std::set<OrbitElementType, CompareType> make_orbit(
    OrbitElementType const &orbit_element, GroupElementIt group_begin,
    GroupElementIt group_end, CompareType compare_f,
    CopyApplyType copy_apply_f) {
  std::set<ElementType, CompareElementType> orbit(compare_f);
  for (; group_begin != group_end; ++group_begin) {
    orbit.emplace(copy_apply_f(*group_begin, orbit_element));
  }
  return orbit;
);

template <typename OrbitElementType, typename GroupElementIt,
          typename CompareType, typename CopyApplyType>
std::vector<std::vector<Index>> make_equivalence_map(
    std::set<OrbitElementType, CompareType> const &orbit,
    GroupElementIt group_begin, GroupElementIt group_end,
    CopyApplyType copy_apply_f) {
  std::vector<std::vector<Index>> equivalence_map;
  equivalence_map.resize(orbit.size());
  Index i = 0;
  for (; group_begin != group_end; ++group_begin) {
    auto it = orbit.find(copy_apply_f(*group_begin, orbit_element));
    if (it == orbit.end()) {
      throw std::runtime_error("Error in make_equivalence_map: failed");
    }
    Index d = std::distance(orbit.begin(), it);
    equivalence_map[d].push_back(i);
    ++i;
  }
  return equivalence_map;
}

}  // namespace group
}  // namespace CASM

#endif
