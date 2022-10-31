#include "casm/configuration/FromStructure.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSet.hh"
#include "casm/configuration/definitions.hh"
#include "casm/configuration/misc.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/StrainConverter.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
namespace config {

FromStructure::FromStructure(std::shared_ptr<Prim const> const &_prim)
    : m_prim(_prim), m_xtal_prim(m_prim->basicstructure) {}

std::runtime_error FromStructure::error(std::string what) const {
  std::string prefix = "Error in FromStructure: ";
  return std::runtime_error(prefix + what);
}

AnisoValTraits FromStructure::get_local_traits_or_throw(std::string key) const {
  try {
    return CASM::config::get_local_traits_or_throw(key);
  } catch (std::exception &e) {
    throw this->error(std::string("Cannot transform local property '") + key +
                      "'");
  }
}

AnisoValTraits FromStructure::get_global_traits_or_throw(
    std::string key) const {
  try {
    return CASM::config::get_local_traits_or_throw(key);
  } catch (std::exception &e) {
    throw this->error(std::string("Cannot transform global property '") + key +
                      "'");
  }
}

Eigen::VectorXd FromStructure::get_Ustrain_vector(
    xtal::SimpleStructure const &mapped_structure) const {
  auto Ustrain_it = mapped_structure.properties.find("Ustrain");
  if (Ustrain_it == mapped_structure.properties.end()) {
    throw this->error("missing required global property 'Ustrain'");
  }
  Eigen::VectorXd Ustrain_vector = Ustrain_it->second;
  if (Ustrain_vector.size() != 6) {
    std::stringstream msg;
    msg << "global property 'Ustrain' size (" << Ustrain_vector.size()
        << ") != 6";
    throw this->error(msg.str());
  }
  return Ustrain_vector;
}

Eigen::Matrix3d FromStructure::get_Ustrain_matrix(
    xtal::SimpleStructure const &mapped_structure) const {
  Eigen::VectorXd Ustrain_vector = get_Ustrain_vector(mapped_structure);
  xtal::StrainConverter Ustrain_converter("Ustrain",
                                          Eigen::MatrixXd::Identity(6, 6));
  return Ustrain_converter.to_E_matrix(Ustrain_vector);
}

std::shared_ptr<Supercell const> FromStructure::make_supercell(
    xtal::SimpleStructure const &mapped_structure) const {
  // --- Construct supercell from lattice and Ustrain ---
  // L_mapped = Ustrain * L_ideal
  Eigen::Matrix3d U = get_Ustrain_matrix(mapped_structure);
  Eigen::Matrix3d const &L_mapped = mapped_structure.lat_column_mat;
  Eigen::Matrix3d L_ideal = U.inverse() * L_mapped;

  double xtal_tol = m_xtal_prim->lattice().tol();
  return std::make_shared<Supercell>(m_prim, xtal::Lattice(L_ideal, xtal_tol));
}

Eigen::MatrixXd const &FromStructure::get_local_property_or_throw(
    std::string key,
    std::map<std::string, Eigen::MatrixXd> const &properties) const {
  auto it = properties.find(key);
  if (it == properties.end()) {
    std::stringstream msg;
    msg << "Missing local property '" << key << "'";
    throw this->error(msg.str());
  }
  return it->second;
}

Eigen::MatrixXd const &FromStructure::get_global_property_or_throw(
    std::string key,
    std::map<std::string, Eigen::MatrixXd> const &properties) const {
  auto it = properties.find(key);
  if (it == properties.end()) {
    std::stringstream msg;
    msg << "Missing global property '" << key << "'";
    throw this->error(msg.str());
  }
  return it->second;
}

void FromStructure::validate_atom_coords_or_throw(
    xtal::SimpleStructure const &mapped_structure,
    std::shared_ptr<Supercell const> const &supercell) const {
  auto const &converter = supercell->unitcellcoord_index_converter;
  auto const &xtal_prim = supercell->prim->basicstructure;
  Index n_sites = converter.total_sites();
  Eigen::MatrixXd R(3, n_sites);
  for (Index l = 0; l < n_sites; ++l) {
    xtal::UnitCellCoord bijk = converter(l);
    R.col(l) = bijk.coordinate(*xtal_prim).const_cart();
  }

  Eigen::MatrixXd disp;
  auto it = mapped_structure.atom_info.properties.find("disp");
  if (it != mapped_structure.atom_info.properties.end()) {
    disp = it->second;
  } else {
    disp = Eigen::MatrixXd::Zero(3, n_sites);
  }

  Eigen::MatrixXd const &coords = mapped_structure.atom_info.coords;
  if (coords.rows() != 3 || coords.cols() != n_sites) {
    std::stringstream msg;
    msg << "atom coords shape (" << coords.rows() << "," << coords.cols()
        << ") does not match the number of supercell sites (" << n_sites << ")";
    throw this->error(msg.str());
  }

  if (!almost_equal(R + disp, coords)) {
    jsonParser json;
    json["supercell_lattice_row_vectors"] =
        supercell->superlattice.superlattice().lat_column_mat().transpose();
    json["supercell_site_coordinates_cart"] = R.transpose();
    json["disp"] = disp.transpose();
    json["displaced_coordinates"] = (R + disp).transpose();
    json["coords"] = coords.transpose();

    fs::path name = error_filename();
    json.write(name);
    std::stringstream msg;
    msg << "atomic displacements ('disp') and atomic coordinates are not "
           "consistent with the supercell site coordinates. See "
        << name << ".";
    throw this->error(msg.str());
  }
}

void FromStructure::validate_atom_names_or_throw(
    std::vector<std::string> const &names, Index n_sites) const {
  if (names.size() != n_sites) {
    std::stringstream msg;
    msg << "the number of atom names (" << names.size()
        << ") does not match the number of sites (" << n_sites << ")";
    throw this->error(msg.str());
  }
}

void FromStructure::validate_local_property_or_throw(
    std::string key, Eigen::MatrixXd const &standard_dof_values,
    Index n_sites) const {
  AnisoValTraits traits = this->get_local_traits_or_throw(key);
  if (standard_dof_values.cols() != n_sites ||
      standard_dof_values.rows() != traits.dim()) {
    std::stringstream msg;
    msg << "atom property '" << key << "' has shape ("
        << standard_dof_values.rows() << "," << standard_dof_values.cols()
        << ") which does not match the expected shape (" << traits.dim() << ","
        << n_sites << ")";
    throw this->error(msg.str());
  }
}

void FromStructure::validate_global_property_or_throw(
    std::string key, Eigen::MatrixXd const &standard_dof_values) const {
  AnisoValTraits traits = this->get_local_traits_or_throw(key);
  if (standard_dof_values.cols() != 1 ||
      standard_dof_values.rows() != traits.dim()) {
    std::stringstream msg;
    msg << "global property '" << key << "' has shape ("
        << standard_dof_values.rows() << "," << standard_dof_values.cols()
        << ") which does not match the expected shape (" << traits.dim()
        << ",1)";
    throw this->error(msg.str());
  }
}

fs::path FromStructure::error_filename() const {
  std::string name = "structure_mapping_error.json";
  if (fs::exists(name)) {
    Index i = 1;
    while (fs::exists(std::string("structure_mapping_error.") +
                      std::to_string(i) + ".json")) {
      ++i;
    }
  }
  return name;
}

/// \brief Constructor
///
/// \param _prim The prim
/// \param _supercells Shared pointer to a SupercellSet. This is an
///     optional method to avoid unnecessarily creating duplicate
///     Supercell. An empty SupercellSet is constructed by default.
FromIsotropicAtomicStructure::FromIsotropicAtomicStructure(
    std::shared_ptr<Prim const> const &_prim,
    std::shared_ptr<SupercellSet> _supercells)
    : FromStructure(_prim), m_supercells(_supercells) {
  if (m_supercells == nullptr) {
    m_supercells = std::make_shared<SupercellSet>(m_prim);
  }
}

/// \brief Construct a configuration with properties from a mapped structure
ConfigurationWithProperties FromIsotropicAtomicStructure::operator()(
    xtal::SimpleStructure const &mapped_structure) {
  // this must be first so other methods can use m_current_supercell
  m_current_supercell =
      m_supercells->insert(this->make_supercell(mapped_structure))
          .first->supercell;

  clexulator::ConfigDoFValues dof_values;
  dof_values.occupation = this->make_occupation(mapped_structure);
  dof_values.local_dof_values = this->make_local_dof_values(mapped_structure);
  dof_values.global_dof_values = this->make_global_dof_values(mapped_structure);

  return ConfigurationWithProperties(
      Configuration(m_current_supercell, dof_values),
      this->make_local_properties(mapped_structure),
      this->make_global_properties(mapped_structure));
}

/// \brief Return shared pointer to all supercells
std::shared_ptr<SupercellSet> FromIsotropicAtomicStructure::supercells() const {
  return m_supercells;
}

std::runtime_error FromIsotropicAtomicStructure::error(std::string what) const {
  std::string prefix = "Error in FromIsotropicAtomicStructure: ";
  return std::runtime_error(prefix + what);
}

Eigen::VectorXi FromIsotropicAtomicStructure::make_occupation(
    xtal::SimpleStructure const &mapped_structure) const {
  auto const &converter = m_current_supercell->unitcellcoord_index_converter;
  Index n_sites = converter.total_sites();
  validate_atom_names_or_throw(mapped_structure.atom_info.names, n_sites);

  Eigen::VectorXi occupation(n_sites);

  // set occupation - by matching atom_info.names to Molecule::name()
  Index l = 0;
  for (std::string const &name : mapped_structure.atom_info.names) {
    xtal::UnitCellCoord bijk = converter(l);
    xtal::Site const &site = m_xtal_prim->basis()[bijk.sublattice()];
    Index s = 0;
    for (auto const &mol : site.occupant_dof()) {
      if (mol.name() == name) {
        break;
      }
    }
    if (s == site.occupant_dof().size()) {
      std::stringstream msg;
      msg << "Failed constructing occupation: atom type '" << name
          << "' is not allowed on site " << l;
      throw this->error(msg.str());
    }
    occupation(l) = s;
    ++l;
  }
  return occupation;
}

std::map<std::string, Eigen::MatrixXd>
FromIsotropicAtomicStructure::make_local_dof_values(
    xtal::SimpleStructure const &mapped_structure) const {
  Index n_sublat = m_xtal_prim->basis().size();
  Index n_unitcells =
      m_current_supercell->unitcell_index_converter.total_sites();
  Index n_sites =
      m_current_supercell->unitcellcoord_index_converter.total_sites();

  std::map<std::string, Eigen::MatrixXd> local_dof_values;
  for (auto const &pair : m_prim->local_dof_info) {
    std::string key = pair.first;
    auto const &dof_info = pair.second;
    Eigen::MatrixXd const &standard_dof_values =
        get_local_property_or_throw(key, mapped_structure.atom_info.properties);
    validate_local_property_or_throw(key, standard_dof_values, n_sites);

    // set DoF, converting from standard to prim basis
    local_dof_values.at(key) = clexulator::local_from_standard_values(
        standard_dof_values, n_sublat, n_unitcells, dof_info);
  }
  return local_dof_values;
}

std::map<std::string, Eigen::MatrixXd>
FromIsotropicAtomicStructure::make_local_properties(
    xtal::SimpleStructure const &mapped_structure) const {
  Index n_sites =
      m_current_supercell->unitcellcoord_index_converter.total_sites();

  std::map<std::string, Eigen::MatrixXd> local_properties;
  for (auto const &atom_property : mapped_structure.atom_info.properties) {
    std::string const &key = atom_property.first;
    Eigen::MatrixXd const &standard_dof_values = atom_property.second;

    // skip DoF
    auto it = m_prim->local_dof_info.find(key);
    if (it == m_prim->local_dof_info.end()) {
      continue;
    }

    // check type and shape
    validate_local_property_or_throw(key, standard_dof_values, n_sites);

    // set as a property
    local_properties.emplace(atom_property);
  }
  return local_properties;
}

std::map<std::string, Eigen::VectorXd>
FromIsotropicAtomicStructure::make_global_dof_values(
    xtal::SimpleStructure const &mapped_structure) const {
  std::map<std::string, Eigen::VectorXd> global_dof_values;
  for (auto const &pair : m_prim->global_dof_info) {
    std::string dof_key = pair.first;
    auto const &dof_info = pair.second;
    bool is_strain = (dof_key.find("strain") != std::string::npos);

    if (is_strain) {
      // get Ustrain and convert to DoF strain metric in prim basis
      xtal::StrainConverter Ustrain_converter("Ustrain",
                                              Eigen::MatrixXd::Identity(6, 6));
      xtal::StrainConverter DoFstrain_converter(
          dof_key, m_prim->global_dof_info.at(dof_key).basis());
      Eigen::VectorXd Ustrain_vector = get_Ustrain_vector(mapped_structure);
      Eigen::Matrix3d F = Ustrain_converter.to_F(Ustrain_vector);
      global_dof_values.at(dof_key) = DoFstrain_converter.from_F(F);
    } else {
      std::string property_key = dof_key;
      Eigen::MatrixXd const &standard_dof_values = get_global_property_or_throw(
          property_key, mapped_structure.properties);
      validate_global_property_or_throw(property_key, standard_dof_values);

      // convert from standard to prim basis
      global_dof_values.at(dof_key) = clexulator::global_from_standard_values(
          standard_dof_values, dof_info);
    }
  }
  return global_dof_values;
}

std::map<std::string, Eigen::VectorXd>
FromIsotropicAtomicStructure::make_global_properties(
    xtal::SimpleStructure const &mapped_structure) const {
  std::map<std::string, Eigen::VectorXd> global_properties;
  for (auto const &global_property : mapped_structure.properties) {
    std::string const &property_key = global_property.first;
    bool is_strain = (property_key.find("strain") != std::string::npos);
    if (is_strain && has_strain_dof(*m_xtal_prim)) {
      continue;
    }
    Eigen::MatrixXd const &standard_dof_values = global_property.second;
    validate_global_property_or_throw(property_key, standard_dof_values);

    // set as a property
    global_properties.emplace(property_key, standard_dof_values);
  }
  return global_properties;
}

}  // namespace config
}  // namespace CASM
