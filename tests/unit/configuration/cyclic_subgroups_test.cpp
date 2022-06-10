#include "casm/configuration/PrimSymInfo.hh"
#include "casm/configuration/basic_symmetry/subgroups.hh"
#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

namespace cyclic_subgroups_test {

Index multiply_f(Index i, Index j) { return (i + j) % 10; }

bool equal_to_f(Index i, Index j) { return i == j; }

std::function<bool(config::basic_symmetry::SubgroupIndices const &)>
make_any_count(
    std::set<config::basic_symmetry::SubgroupOrbit> const &cyclic_subgroups) {
  return [&](config::basic_symmetry::SubgroupIndices const &subgroup) {
    for (auto const &orbit : cyclic_subgroups) {
      if (orbit.count(subgroup)) {
        return true;
      }
    }
    return false;
  };
}

}  // namespace cyclic_subgroups_test

TEST(CyclicSubgroupsTest, Test1) {
  using namespace config::basic_symmetry;
  using namespace cyclic_subgroups_test;
  std::vector<Index> elements;
  for (Index i = 0; i < 10; ++i) {
    elements.push_back(i);
  }
  Group<Index> group = make_group(elements, multiply_f, equal_to_f);
  std::set<SubgroupOrbit> cyclic_subgroups = make_cyclic_subgroups(group);

  EXPECT_EQ(cyclic_subgroups.size(), 4);
  auto any_count = make_any_count(cyclic_subgroups);
  EXPECT_EQ(any_count({0}), 1);
  EXPECT_EQ(any_count({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 1);
  EXPECT_EQ(any_count({0, 2, 4, 6, 8}), 1);
  EXPECT_EQ(any_count({0, 5}), 1);
}

TEST(CyclicSubgroupsTest, Test2) {
  using namespace config::basic_symmetry;
  using namespace cyclic_subgroups_test;
  config::PrimSymInfo prim_sym_info(test::FCC_binary_prim());
  std::set<SubgroupOrbit> cyclic_subgroups =
      make_cyclic_subgroups(*prim_sym_info.factor_group);

  EXPECT_EQ(cyclic_subgroups.size(), 10);
  auto any_count = make_any_count(cyclic_subgroups);
  EXPECT_EQ(any_count({0}), 1);
  EXPECT_EQ(any_count({0, 1, 4, 17}), 1);
  EXPECT_EQ(any_count({0, 7, 11}), 1);
  EXPECT_EQ(any_count({0, 7, 11, 33, 37, 47}), 1);
  EXPECT_EQ(any_count({0, 15}), 1);
  EXPECT_EQ(any_count({0, 17}), 1);
  EXPECT_EQ(any_count({0, 17, 41, 44}), 1);
  EXPECT_EQ(any_count({0, 24}), 1);
  EXPECT_EQ(any_count({0, 26}), 1);
  EXPECT_EQ(any_count({0, 47}), 1);
}

TEST(AllSubgroupsTest, Test1) {
  using namespace config::basic_symmetry;
  using namespace cyclic_subgroups_test;
  std::vector<Index> elements;
  for (Index i = 0; i < 10; ++i) {
    elements.push_back(i);
  }
  Group<Index> group = make_group(elements, multiply_f, equal_to_f);
  std::set<SubgroupOrbit> all_subgroups = make_all_subgroups(group);

  EXPECT_EQ(all_subgroups.size(), 4);
  auto any_count = make_any_count(all_subgroups);
  EXPECT_EQ(any_count({0}), 1);
  EXPECT_EQ(any_count({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}), 1);
  EXPECT_EQ(any_count({0, 2, 4, 6, 8}), 1);
  EXPECT_EQ(any_count({0, 5}), 1);
}

TEST(AllSubgroupsTest, Test2) {
  using namespace config::basic_symmetry;
  using namespace cyclic_subgroups_test;
  config::PrimSymInfo prim_sym_info(test::FCC_binary_prim());
  std::set<SubgroupOrbit> all_subgroups =
      make_all_subgroups(*prim_sym_info.factor_group);

  EXPECT_EQ(all_subgroups.size(), 33);
  auto any_count = make_any_count(all_subgroups);
  EXPECT_EQ(any_count({
                0,
            }),
            1);
  EXPECT_EQ(any_count({
                0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            }),
            1);
  EXPECT_EQ(any_count({
                0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                1,
                4,
                17,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                1,
                4,
                17,
                18,
                21,
                22,
                23,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                1,
                4,
                17,
                18,
                21,
                22,
                23,
                26,
                27,
                30,
                31,
                32,
                41,
                44,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                1,
                4,
                17,
                26,
                41,
                44,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                1,
                4,
                17,
                27,
                30,
                31,
                32,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                17,
                21,
                22,
            }),
            1);
  EXPECT_EQ(any_count({
                0,  7,  8,  9,  10, 11, 12, 13, 14, 17, 21, 22,
                24, 25, 27, 28, 29, 32, 41, 42, 43, 44, 45, 46,
            }),
            1);
  EXPECT_EQ(any_count({
                0,  7,  8,  9,  10, 11, 12, 13, 14, 17, 21, 22,
                26, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                11,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                11,
                18,
                19,
                20,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                11,
                18,
                19,
                20,
                27,
                28,
                29,
                33,
                37,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                11,
                27,
                28,
                29,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                7,
                11,
                33,
                37,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
                20,
                21,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
                20,
                21,
                24,
                29,
                30,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
                20,
                21,
                26,
                31,
                42,
                45,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
                24,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                15,
                29,
                30,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                21,
                22,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                21,
                22,
                24,
                29,
                42,
                45,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                21,
                22,
                26,
                30,
                31,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                26,
                47,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                27,
                32,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                30,
                31,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                17,
                41,
                44,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                24,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                26,
            }),
            1);
  EXPECT_EQ(any_count({
                0,
                47,
            }),
            1);
}
