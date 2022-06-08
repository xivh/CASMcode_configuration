#include "casm/configuration/DoFSpaceAxisInfo.hh"

#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/LinearIndexConverter.hh"

namespace CASM {
namespace config {

namespace DoFSpace_impl {

void throw_if_missing_local_dof_requirements(
    DoFKey const &dof_key,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  if (!transformation_matrix_to_super.has_value() || !sites.has_value()) {
    std::stringstream msg;
    msg << "Error: local DoF '" << dof_key
        << "' require transformation_matrix_to_super and sites" << std::endl;
    throw std::runtime_error(msg.str());
  }
}

}  // namespace DoFSpace_impl

/// \brief Make DoFSpace axis glossary, axis site index, and axis dof component
DoFSpaceAxisInfo::DoFSpaceAxisInfo(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  if (AnisoValTraits(dof_key).global()) {
    glossary.clear();
    site_index = std::nullopt;
    dof_component = std::nullopt;
    basis_row_index = std::nullopt;

    // Global DoF, glossary comes straight from the DoF
    glossary = component_descriptions(prim.global_dof(dof_key));
  } else {
    using namespace DoFSpace_impl;
    throw_if_missing_local_dof_requirements(
        dof_key, transformation_matrix_to_super, sites);

    glossary.clear();
    site_index = std::vector<Index>();
    dof_component = std::vector<Index>();
    basis_row_index = std::vector<std::vector<Index>>();

    xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
        *transformation_matrix_to_super, prim.basis().size());
    basis_row_index->resize(unitcellcoord_index_converter.total_sites());

    Index _basis_row_index = 0;
    // Generate full glossary for all active sites of the config_region
    for (Index l : *sites) {
      Index sublattice_index = unitcellcoord_index_converter(l).sublattice();
      xtal::Site const &site = prim.basis()[sublattice_index];
      if (dof_key == "occ") {
        Index i = 0;
        for (auto const &molecule : site.occupant_dof()) {
          glossary.push_back("occ[" + std::to_string(l + 1) + "][" +
                             molecule.name() + "]");
          dof_component->push_back(i);
          site_index->push_back(l);
          (*basis_row_index)[l].push_back(_basis_row_index);
          ++i;
          ++_basis_row_index;
        }
      } else if (site.has_dof(dof_key)) {
        std::vector<std::string> tdescs =
            component_descriptions(site.dof(dof_key));
        Index i = 0;
        for (std::string const &desc : tdescs) {
          glossary.push_back(desc + "[" + std::to_string(l + 1) + "]");
          dof_component->push_back(i);
          site_index->push_back(l);
          (*basis_row_index)[l].push_back(_basis_row_index);
          ++i;
          ++_basis_row_index;
        }
      }
      if (!site.has_dof(dof_key)) continue;
    }
  }
}

}  // namespace config
}  // namespace CASM
