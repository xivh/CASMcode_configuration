#ifndef CASM_config_Prim
#define CASM_config_Prim

#include "casm/configuration/PrimSymInfo.hh"
#include "casm/configuration/definitions.hh"
#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace config {

/// \brief Species the primitive crystal structure (lattice and basis) and
/// allowed degrees of freedom (DoF), and also symmetry representations
/// used for all configurations with the same prim. All members are const.
struct Prim {
  Prim(BasicStructure const &_basicstructure);

  /// \brief The BasicStructure specifies the primitive crystal structure
  /// (lattice and basis) and allowed degrees of freedom (DoF)
  BasicStructure const basicstructure;

  /// \brief Global DoFSet information, defines bases for global DoF
  std::map<DoFKey, xtal::DoFSet> const global_dof_info;

  /// \brief Local DoFSet information, defines bases on each sublattice for
  ///     local DoF
  std::map<DoFKey, std::vector<xtal::SiteDoFSet>> const local_dof_info;

  /// \brief Holds symmetry representations used for all configurations with
  /// the same prim structure
  PrimSymInfo const sym_info;
};

}  // namespace config
}  // namespace CASM

#endif
