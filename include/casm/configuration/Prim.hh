#ifndef CASM_config_Prim
#define CASM_config_Prim

#include "casm/configuration/PrimMagspinInfo.hh"
#include "casm/configuration/PrimSymInfo.hh"
#include "casm/configuration/definitions.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace config {

/// \brief Species the primitive crystal structure (lattice and basis) and
/// allowed degrees of freedom (DoF), and also symmetry representations
/// used for all configurations with the same prim. All members are const.
struct Prim {
  /// \brief Constructor
  Prim(std::shared_ptr<BasicStructure const> const &_basicstructure);

  /// \brief Construct using factor group in given order
  Prim(std::vector<xtal::SymOp> const &factor_group_elements,
       std::shared_ptr<BasicStructure const> const &_basicstructure);

  /// \brief The BasicStructure specifies the primitive crystal structure
  /// (lattice and basis) and allowed degrees of freedom (DoF)
  std::shared_ptr<BasicStructure const> const basicstructure;

  /// \brief Global DoFSet information, defines bases for global DoF
  std::map<DoFKey, xtal::DoFSet> const global_dof_info;

  /// \brief Local DoFSet information, defines bases on each sublattice for
  ///     local DoF
  std::map<DoFKey, std::vector<xtal::SiteDoFSet>> const local_dof_info;

  /// \brief Checks that the prim allows 1 or more occupants on each site,
  ///     all occupants have a single atom, and there are no Molecule
  ///     properties, only AtomPosition properties
  bool is_atomic;

  /// \brief Holds symmetry representations used for all configurations with
  /// the same prim structure
  PrimSymInfo const sym_info;

  /// \brief Holds information about allowed occupants with discrete magspin or
  ///     continuous magspin DoF
  PrimMagspinInfo const magspin_info;
};

inline std::shared_ptr<Prim const> make_shared_prim(
    xtal::BasicStructure const &basicstructure) {
  return std::make_shared<config::Prim const>(
      std::make_shared<xtal::BasicStructure const>(basicstructure));
}

inline std::shared_ptr<Prim const> make_shared_prim(
    std::shared_ptr<xtal::BasicStructure const> const &basicstructure) {
  return std::make_shared<config::Prim const>(basicstructure);
}

}  // namespace config
}  // namespace CASM

#endif
