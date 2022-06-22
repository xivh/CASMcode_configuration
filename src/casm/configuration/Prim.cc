#include "casm/configuration/Prim.hh"

#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"

namespace CASM {
namespace config {

Prim::Prim(std::shared_ptr<BasicStructure const> const &_basicstructure)
    : basicstructure(_basicstructure),
      global_dof_info(clexulator::make_global_dof_info(*basicstructure)),
      local_dof_info(clexulator::make_local_dof_info(*basicstructure)),
      sym_info(*basicstructure) {}

}  // namespace config
}  // namespace CASM
