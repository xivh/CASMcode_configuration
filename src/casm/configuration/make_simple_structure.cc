#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/StrainConverter.hh"

namespace CASM {
namespace config {

/// \brief Convert a Configuration to a SimpleStructure
///
/// \param configuration The configuration being converted to
///     a SimpleStructure. If "disp" is a DoF, it is included
///     in SimpleStructure::atom_info.coords, but not included
///     in SimpleStructure::atom_info.properties. If strain is
///     a DoF, it is included in SimpleStructure::atom_info.coords
///     and SimpleStructure::lat_column_mat, but not included in
///     SimpleStructure::properties. Other DoF are copied to
///     SimpleStructure::atom_info.properties or
///     SimpleStructure::properties.
/// \param local_properties If provided, these are copied to
///     SimpleStructure::atom_info.properties. If "disp" is included,
///     it is included in SimpleStructure::atom_info.coords.
/// \param global_properties If provided, these are copied to
///     SimpleStructure::properties. If "Ustrain" is included,
///     it is included in SimpleStructure::lat_column_mat.
///
/// \returns simple_structure A SimpleStructure representation
///     of the configuration
///
/// Notes:
/// - This is a simpler method which accepts configurations from
///   prim with only atomic occupants
/// - Deformation gradient applied to ideal lattice vectors is obtained from:
///   1) a strain DoF (this is not copied to structure.properties)
///   2) "Ustrain" in global_properties (this is copied to structure.properties)
///   3) the identify matrix, I
/// -
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration,
    std::map<std::string, Eigen::MatrixXd> const &local_properties,
    std::map<std::string, Eigen::VectorXd> const &global_properties) {
  // references
  auto const &supercell = *configuration.supercell;
  auto const &prim = *supercell.prim;
  auto const &converter = supercell.unitcellcoord_index_converter;
  Index n_sites = converter.total_sites();
  Index N_sublat = prim.basicstructure->basis().size();
  Index N_unitcells = supercell.unitcell_index_converter.total_sites();
  auto const &dof_values = configuration.dof_values;
  auto const &occupation = dof_values.occupation;
  auto const &local_dof_values = dof_values.local_dof_values;
  auto const &global_dof_values = dof_values.global_dof_values;

  // validate no multi-atom molecules
  std::vector<xtal::Molecule> molecule_list =
      xtal::struc_molecule(*prim.basicstructure);
  for (auto const &mol : molecule_list) {
    if (mol.size() > 1) {
      throw std::runtime_error(
          "Error in make_simple_structure:  multi-atom molecules are not "
          "allowed");
    }
  }

  // set m_ideal_lat_column_mat
  Eigen::Matrix3d ideal_lat_column_mat =
      supercell.superlattice.superlattice().lat_column_mat();

  // get ideal site coordinates
  Eigen::MatrixXd coords(3, n_sites);
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    coords.col(l) = bijk.coordinate(*prim.basicstructure).const_cart();
  }

  // get atom type names
  std::vector<std::string> names;
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    Index b = converter(l).sublattice();
    auto const &site = prim.basicstructure->basis()[b];
    names.push_back(site.occupant_dof()[occupation(l)].name());
  }

  // get deformation gradient, F
  Eigen::Matrix3d F;
  DoFKey strain_dof_key;
  DoFKey strain_property_key;
  if (has_strain_dof(*prim.basicstructure)) {
    strain_dof_key = get_strain_dof_key(*prim.basicstructure);
    xtal::StrainConverter DoFstrain_converter(
        strain_dof_key, prim.global_dof_info.at(strain_dof_key).basis());
    F = DoFstrain_converter.to_F(global_dof_values.at(strain_dof_key));
  } else if (global_properties.count("Ustrain")) {
    strain_property_key = "Ustrain";
    xtal::StrainConverter Ustrain_converter("Ustrain",
                                            Eigen::MatrixXd::Identity(6, 6));
    F = Ustrain_converter.to_F(global_properties.at(strain_property_key));
  } else {
    F = Eigen::Matrix3d::Identity();
  }

  // get displacements
  Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(3, n_sites);
  if (local_dof_values.count("disp")) {
    disp = clexulator::local_to_standard_values(local_dof_values.at("disp"),
                                                N_sublat, N_unitcells,
                                                prim.local_dof_info.at("disp"));
  } else if (local_properties.count("disp")) {
    disp = local_properties.at("disp");
  }

  // construct SimpleStructure
  xtal::SimpleStructure structure;

  structure.lat_column_mat = ideal_lat_column_mat;
  structure.atom_info.resize(names.size());
  structure.atom_info.names = names;
  for (Index i = 0; i < names.size(); ++i) {
    structure.atom_info.coords = coords + disp;
  }

  // copy local_dof (excluding disp)
  for (auto const &local_dof : local_dof_values) {
    if (local_dof.first != "disp") {
      structure.atom_info.properties.insert(local_dof);
    }
  }

  // copy global_dof (excluding strain)
  for (auto const &global_dof : global_dof_values) {
    if (global_dof.first != strain_dof_key) {
      structure.properties.insert(global_dof);
    }
  }

  // copy local_properties
  for (auto const &local_property : local_properties) {
    structure.atom_info.properties.insert(local_property);
  }

  // copy global_properties
  for (auto const &global_property : global_properties) {
    structure.properties.insert(global_property);
  }

  structure.deform_coords(F);
  return structure;
}

}  // namespace config
}  // namespace CASM
