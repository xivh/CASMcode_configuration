// The `casm/configuration` module enables apply symmetry to ConfigDoFValues.
//
// Primarily, this purpose of this module is to provide:
//
// - `struct Configuration`: Data structure encapsulating configuration degree
// of freedom (DoF) values and all information necessary to apply symmetry
// operations. This data structure contains:
//   - `ConfigDoFValues dof_values`: The raw DoF values
//   - `std::shared_ptr<Supercell const> supercell`: Specifies all the
//   structural and symmetry information common for all configurations with the
//   same supercell. All members are const.
// - `class SupercellOpIterator`: Class that enables iterating over all symmetry
// operations that are compatible with a particular supercell and applying them
// to a Configuration.
// - Methods for applying symmetry operations:
//   - `ConfigDoFValues &apply(
//          SupercellOpIterator const &op,
//          ConfigDoFValues &dof_values)`;
//   - `ConfigDoFValues copy_apply(
//          SupercellOpIterator const &op,
//          ConfigDoFValues dof_values);`
// - Methods for comparing configurations and finding canonical forms:
//  - `class ConfigCompare`: Provides "less than" comparison of Configuration
//  - `class ConfigIsEquivalent`: Provides "equal to" comparison of
//  Configuration
//  - `is_canonical`, `canonical_form`, `to_canonical`, `from_canonical`:
//     Methods that provide checking and finding Configuration canonical forms
//     (the "greatest" of all symmetrically equivalent configuration).
//
// Allowed dependencies:
// - CASMcode_global
// - CASMcode_crystallography
// - CASMcode_clexulator
// - CASMcode_configuration/group
// - CASMcode_configuration/sym_info
//
// Not allowed:
// - CASMcode_configuration/clusterography
// - CASMcode_configuration/occ_events

#ifndef CASM_config_definitions
#define CASM_config_definitions

#include <memory>
#include <optional>
#include <set>
#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace clexulator {
struct ConfigDoFValues;
struct DoFSpace;
}  // namespace clexulator

namespace xtal {
class BasicStructure;
class Lattice;
class SimpleStructure;
class Superlattice;
struct SymOp;
class UnitCell;
class UnitCellCoord;
struct UnitCellCoordRep;
class UnitCellIndexConverter;
class UnitCellCoordIndexConverter;
}  // namespace xtal

namespace group {
template <typename ElementType>
struct Group;
}

namespace config {

using clexulator::ConfigDoFValues;
using xtal::BasicStructure;
using xtal::Lattice;
using xtal::SimpleStructure;
using xtal::Superlattice;
using xtal::SymOp;
using xtal::UnitCell;
using xtal::UnitCellCoord;
using xtal::UnitCellCoordRep;

struct Configuration;
struct Prim;
struct PrimSymInfo;
struct Supercell;
struct SupercellSymInfo;
class SupercellSymOp;

typedef long Index;
typedef std::string DoFKey;

/// \brief A group::Group of xtal::SymOp
typedef group::Group<SymOp> SymGroup;

template <typename T>
T const &throw_if_equal_to_nullptr(T const &t, std::string message);

// --- Inline definitions ---

template <typename T>
T const &throw_if_equal_to_nullptr(T const &t, std::string message) {
  if (t == nullptr) {
    throw std::runtime_error(message);
  }
  return t;
}

}  // namespace config
}  // namespace CASM

#endif
