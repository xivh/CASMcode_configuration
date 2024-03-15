#include "casm/configuration/make_simple_structure.hh"

#include <vector>  // see https://github.com/prisms-center/CASMcode_clexulator/issues/19

#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/StrainConverter.hh"

namespace CASM {
namespace config {

/// \brief (deprecated) Convert a Configuration to a SimpleStructure
///
/// This method is deprecated in favor of ToAtomicStructure.
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
/// \param atom_type_naming_method Specifies how to set atom_info.names.
/// \param excluded_species Specifies names that should not be
///     included in the resulting xtal::SimpleStructure.
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
/// - Options for `atom_type_naming_method` are:
///   - "chemical_name" (default): use ``mol.name()``, where `mol`
///     is the `xtal::Molecule` occupying the site
///   - "orientation_name": use ``unique_names[b][s]``, where `unique_names` is
///     the unique occupant names obtained from
///     ``xtal::BasicStructure::unique_names``, `b` is the sublattice index of
///     the site, and `s` is the occupation index.
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration,
    std::map<std::string, Eigen::MatrixXd> const &local_properties,
    std::map<std::string, Eigen::VectorXd> const &global_properties,
    std::string atom_type_naming_method,
    std::set<std::string> excluded_species) {
  // references
  auto const &supercell = *configuration.supercell;
  auto const &prim = *supercell.prim;
  auto const &basis = prim.basicstructure->basis();
  auto const &converter = supercell.unitcellcoord_index_converter;
  Index n_sites = converter.total_sites();
  Index N_sublat = basis.size();
  Index N_unitcells = supercell.unitcell_index_converter.total_sites();
  auto const &dof_values = configuration.dof_values;
  auto const &occupation = dof_values.occupation;
  auto const &local_dof_values = dof_values.local_dof_values;
  auto const &global_dof_values = dof_values.global_dof_values;

  // validate no multi-atom molecules
  if (!prim.is_atomic) {
    throw std::runtime_error(
        "Error in make_simple_structure: not an atomic structure");
  }

  // validate only continuous or discrete magspin
  if (prim.magspin_info.has_discrete_atomic_magspin_occupants &&
      prim.magspin_info.has_continuous_magspin_dof) {
    throw std::runtime_error(
        "Error in make_simple_structure: prim has continuous and discrete "
        "magspin");
  }

  // validate occupation indices
  for (Index l = 0; l < occupation.size(); ++l) {
    auto const &site = basis[converter(l).sublattice()];
    int s = occupation(l);
    if (s < 0 || s >= site.occupant_dof().size()) {
      std::stringstream msg;
      msg << "Error in make_simple_structure: invalid occupation=" << s
          << " at linear_site_index=" << l << ".";
      throw std::runtime_error(msg.str());
    }
  }

  // set m_ideal_lat_column_mat
  Eigen::Matrix3d ideal_lat_column_mat =
      supercell.superlattice.superlattice().lat_column_mat();

  // get atom type names
  std::vector<std::string> names;
  std::vector<Index> site_index;
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    Index b = converter(l).sublattice();
    int s = occupation(l);
    std::string name;
    if (atom_type_naming_method == "orientation_name") {
      name = prim.basicstructure->unique_names()[b][s];
    } else if (atom_type_naming_method == "chemical_name") {
      name = basis[b].occupant_dof()[s].name();
    } else {
      std::stringstream msg;
      msg << "Error in make_simple_structure: invalid atom_type_naming_method='"
          << atom_type_naming_method << "'";
      throw std::runtime_error(msg.str());
    }
    if (!excluded_species.count(name)) {
      names.push_back(name);
      site_index.push_back(l);
    }
  }

  // get ideal site coordinates (all sites)
  Eigen::MatrixXd coords(3, n_sites);
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    coords.col(l) = bijk.coordinate(*prim.basicstructure).const_cart();
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

  // get displacements (all sites)
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
  structure.atom_info.coords = Eigen::MatrixXd::Zero(3, site_index.size());
  for (Index i = 0; i < site_index.size(); ++i) {
    Index l = site_index[i];
    structure.atom_info.coords.col(i) = coords.col(l) + disp.col(l);
  }

  // copy local_dof (excluding disp)
  for (auto const &local_dof : local_dof_values) {
    std::string key = local_dof.first;
    Eigen::MatrixXd const &value = local_dof.second;
    Eigen::MatrixXd new_value =
        Eigen::MatrixXd::Zero(value.rows(), site_index.size());
    for (Index i = 0; i < site_index.size(); ++i) {
      Index l = site_index[i];
      new_value.col(i) = value.col(l);
    }
    if (local_dof.first != "disp") {
      structure.atom_info.properties.emplace(key, new_value);
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
    std::string key = local_property.first;
    Eigen::MatrixXd const &value = local_property.second;
    Eigen::MatrixXd new_value =
        Eigen::MatrixXd::Zero(value.rows(), site_index.size());
    for (Index i = 0; i < site_index.size(); ++i) {
      Index l = site_index[i];
      new_value.col(i) = value.col(l);
    }
    structure.atom_info.properties.emplace(key, new_value);
  }

  // copy global_properties
  for (auto const &global_property : global_properties) {
    structure.properties.insert(global_property);
  }

  // selectivedynamics values as atom_info.properties.at("selectivedynamics")
  bool has_selectivedynamics = false;
  for (auto const &site : basis) {
    for (auto const &mol : site.occupant_dof()) {
      auto const &properties = mol.atom(0).properties();
      if (properties.find("selectivedynamics") != properties.end()) {
        has_selectivedynamics = true;
      }
    }
  }
  if (has_selectivedynamics) {
    std::string key = "selectivedynamics";
    AnisoValTraits traits(key);
    Eigen::VectorXd default_atom_value = Eigen::VectorXd::Zero(traits.dim());
    Eigen::MatrixXd structure_values =
        Eigen::MatrixXd::Zero(traits.dim(), names.size());
    for (Index i = 0; i < site_index.size(); ++i) {
      Index l = site_index[i];
      Index b = converter(l).sublattice();
      auto const &site = prim.basicstructure->basis()[b];
      auto const &mol = site.occupant_dof()[occupation(l)];
      auto const &properties = mol.atom(0).properties();
      auto property_it = properties.find(key);
      Eigen::VectorXd atom_value = default_atom_value;
      if (property_it != properties.end()) {
        atom_value = property_it->second.value();
      }
      structure_values.col(i) = atom_value;
    }
    structure.atom_info.properties.emplace(key, structure_values);
  }

  // discrete magspin values as atom_info.properties.at(magspin_key)
  if (prim.magspin_info.has_discrete_atomic_magspin_occupants) {
    if (!prim.magspin_info.discrete_atomic_magspin_key.has_value()) {
      throw std::runtime_error(
          "Error in make_simple_structure: "
          "prim.magspin_info.discrete_atomic_magspin_key has no value");
    }
    std::string key = prim.magspin_info.discrete_atomic_magspin_key.value();
    AnisoValTraits traits(key);
    Eigen::VectorXd default_atom_value = Eigen::VectorXd::Zero(traits.dim());
    Eigen::MatrixXd structure_values =
        Eigen::MatrixXd::Zero(traits.dim(), names.size());
    for (Index i = 0; i < site_index.size(); ++i) {
      Index l = site_index[i];
      Index b = converter(l).sublattice();
      auto const &site = prim.basicstructure->basis()[b];
      auto const &mol = site.occupant_dof()[occupation(l)];
      auto const &properties = mol.atom(0).properties();
      auto property_it = properties.find(key);
      Eigen::VectorXd atom_value = default_atom_value;
      if (property_it != properties.end()) {
        atom_value = property_it->second.value();
      }
      structure_values.col(i) = atom_value;
    }
    structure.atom_info.properties.emplace(key, structure_values);
  }

  structure.deform_coords(F);
  return structure;
}

/// \brief Constructor
///
/// \param atom_type_naming_method Specifies how to set atom_info.names.
/// \param excluded_species Specifies names that should not be
///     included in the resulting xtal::SimpleStructure.
///
/// - Options for `atom_type_naming_method` are:
///   - "chemical_name" (default): use ``mol.name()``, where `mol`
///     is the `xtal::Molecule` occupying the site
///   - "orientation_name": use ``unique_names[b][s]``, where `unique_names` is
///     the unique occupant names obtained from
///     ``xtal::BasicStructure::unique_names``, `b` is the sublattice index of
///     the site, and `s` is the occupation index.
ToAtomicStructure::ToAtomicStructure(std::string atom_type_naming_method,
                                     std::set<std::string> excluded_species)
    : m_atom_type_naming_method(atom_type_naming_method),
      m_excluded_species(excluded_species) {}

/// \brief Convert a ConfigurationWithProperties to a SimpleStructure
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
/// \returns simple_structure A SimpleStructure representation
///     of the configuration
xtal::SimpleStructure ToAtomicStructure::operator()(
    ConfigurationWithProperties const &configuration_with_properties) {
  auto const &x = configuration_with_properties;
  return make_simple_structure(x.configuration, x.local_properties,
                               x.global_properties, m_atom_type_naming_method,
                               m_excluded_species);
}

/// \brief Convert a Configuration to a SimpleStructure
xtal::SimpleStructure ToAtomicStructure::operator()(
    Configuration const &configuration,
    std::map<std::string, Eigen::MatrixXd> const &local_properties,
    std::map<std::string, Eigen::VectorXd> const &global_properties) {
  return make_simple_structure(configuration, local_properties,
                               global_properties, m_atom_type_naming_method,
                               m_excluded_species);
}

}  // namespace config
}  // namespace CASM
