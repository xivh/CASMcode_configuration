#include "casm/configuration/PrimSymInfo.hh"

#include "gtest/gtest.h"
#include "teststructures.hh"

using namespace CASM;

TEST(PrimSymInfoTest, Test1) {
  PrimSymInfo prim_sym_info(test::FCC_binary_prim());

  CHECK_EQ(prim_sym_info.factor_group.size(), 48);
}
