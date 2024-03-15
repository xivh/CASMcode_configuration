#include "casm/configuration/Prim.hh"

#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace config {

namespace {

bool _is_atomic(BasicStructure const &structure) {
  /// Must have:
  /// - 1 or more occupants on each site
  /// - Each occupant has 1 atom and no molecule properties
  /// - Each atom coordinate must be [0., 0., 0.]
  for (auto const &site : structure.basis()) {
    if (site.occupant_dof().size() == 0) {
      return false;
    }
    double xtal_tol = structure.lattice().tol();
    Eigen::Vector3d zero = Eigen::Vector3d::Zero();
    for (auto const &mol : site.occupant_dof()) {
      if (mol.size() != 1 || !mol.properties().empty()) {
        return false;
      }
      if (!almost_equal(mol.atom(0).cart(), zero, xtal_tol)) {
        return false;
      }
    }
  }
  return true;
}

void _validate_unique_names(BasicStructure const &structure) {
  if (structure.unique_names().size() != structure.basis().size()) {
    throw std::runtime_error(
        "Error in config::Prim constructor: invalid unique_names");
  }
  for (Index b = 0; b < structure.basis().size(); ++b) {
    if (structure.unique_names()[b].size() !=
        structure.basis()[b].occupant_dof().size()) {
      throw std::runtime_error(
          "Error in config::Prim constructor: invalid unique_names");
    }
  }
}

}  // namespace

/// \brief Constructor
Prim::Prim(std::shared_ptr<BasicStructure const> const &_basicstructure)
    : basicstructure(throw_if_equal_to_nullptr(
          _basicstructure,
          "Error in Prim constructor: _basicstructure == nullptr")),
      global_dof_info(clexulator::make_global_dof_info(*basicstructure)),
      local_dof_info(clexulator::make_local_dof_info(*basicstructure)),
      is_atomic(_is_atomic(*basicstructure)),
      sym_info(*basicstructure),
      magspin_info(*basicstructure) {
  _validate_unique_names(*basicstructure);
}

/// \brief Construct using factor group in given order
Prim::Prim(std::vector<xtal::SymOp> const &factor_group_elements,
           std::shared_ptr<BasicStructure const> const &_basicstructure)
    : basicstructure(throw_if_equal_to_nullptr(
          _basicstructure,
          "Error in Prim constructor: _basicstructure == nullptr")),
      global_dof_info(clexulator::make_global_dof_info(*basicstructure)),
      local_dof_info(clexulator::make_local_dof_info(*basicstructure)),
      is_atomic(_is_atomic(*basicstructure)),
      sym_info(factor_group_elements, *basicstructure),
      magspin_info(*basicstructure) {
  _validate_unique_names(*basicstructure);
}

}  // namespace config
}  // namespace CASM
