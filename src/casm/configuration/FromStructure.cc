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
  auto validate_f = [this](Eigen::VectorXd const &Xstrain_vector,
                           std::string property_key) {
    if (Xstrain_vector.size() != 6) {
      std::stringstream msg;
      msg << "global property \"" << property_key << "\" size ("
          << Xstrain_vector.size() << ") != 6";
      throw this->error(msg.str());
    }
  };

  auto Ustrain_it = mapped_structure.properties.find("Ustrain");
  if (Ustrain_it == mapped_structure.properties.end()) {
    // throw this->error("missing required global property 'Ustrain'");

    for (auto const &global_property : mapped_structure.properties) {
      std::string const &property_key = global_property.first;
      bool is_strain = (property_key.find("strain") != std::string::npos);
      if (!is_strain) {
        continue;
      }

      // Get mapped_structure strain property and convert to Ustrain
      xtal::StrainConverter Xstrain_converter(property_key,
                                              Eigen::MatrixXd::Identity(6, 6));
      xtal::StrainConverter Ustrain_converter("Ustrain",
                                              Eigen::MatrixXd::Identity(6, 6));
      Eigen::VectorXd Xstrain_vector = global_property.second;
      validate_f(Xstrain_vector, property_key);
      Eigen::Matrix3d F = Xstrain_converter.to_F(Xstrain_vector);
      return Ustrain_converter.from_F(F);
    }

    // No strain found
    Eigen::VectorXd Ustrain_vector(6);
    Ustrain_vector << 1., 1., 1., 0., 0., 0.;
    return Ustrain_vector;
  }
  Eigen::VectorXd Ustrain_vector = Ustrain_it->second;
  validate_f(Ustrain_vector, "Ustrain");
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

std::map<std::string, Eigen::MatrixXd>
FromStructure::default_make_local_dof_values(
    xtal::SimpleStructure const &mapped_structure,
    std::shared_ptr<Supercell const> current_supercell) const {
  Index n_sublat = m_xtal_prim->basis().size();
  Index n_unitcells = current_supercell->unitcell_index_converter.total_sites();
  Index n_sites =
      current_supercell->unitcellcoord_index_converter.total_sites();

  std::map<std::string, Eigen::MatrixXd> local_dof_values;
  for (auto const &pair : m_prim->local_dof_info) {
    std::string key = pair.first;
    auto const &dof_info = pair.second;
    Eigen::MatrixXd const &standard_dof_values =
        get_local_property_or_throw(key, mapped_structure.atom_info.properties);
    validate_local_property_or_throw(key, standard_dof_values, n_sites);

    // set DoF, converting from standard to prim basis
    local_dof_values[key] = clexulator::local_from_standard_values(
        standard_dof_values, n_sublat, n_unitcells, dof_info);
  }
  return local_dof_values;
}

std::map<std::string, Eigen::MatrixXd>
FromStructure::default_make_local_properties(
    xtal::SimpleStructure const &mapped_structure,
    std::shared_ptr<Supercell const> current_supercell,
    std::set<std::string> excluded) const {
  Index n_sites =
      current_supercell->unitcellcoord_index_converter.total_sites();

  std::map<std::string, Eigen::MatrixXd> local_properties;
  for (auto const &atom_property : mapped_structure.atom_info.properties) {
    std::string const &key = atom_property.first;
    Eigen::MatrixXd const &standard_dof_values = atom_property.second;

    // skip excluded
    if (excluded.count(key)) {
      continue;
    }

    // skip DoF
    auto it = m_prim->local_dof_info.find(key);
    if (it != m_prim->local_dof_info.end()) {
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
FromStructure::default_make_global_dof_values(
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
      global_dof_values[dof_key] = DoFstrain_converter.from_F(F);
    } else {
      std::string property_key = dof_key;
      Eigen::MatrixXd const &standard_dof_values = get_global_property_or_throw(
          property_key, mapped_structure.properties);
      validate_global_property_or_throw(property_key, standard_dof_values);

      // convert from standard to prim basis
      global_dof_values[dof_key] = clexulator::global_from_standard_values(
          standard_dof_values, dof_info);
    }
  }
  return global_dof_values;
}

std::map<std::string, Eigen::VectorXd>
FromStructure::default_make_global_properties(
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

  // collect selectivedynamics info
  std::string selectivedynamics_key = "selectivedynamics";
  bool has_selectivedynamics =
      mapped_structure.atom_info.properties.count(selectivedynamics_key);
  Eigen::VectorXd default_selectivedynamics_atom_value;
  Eigen::MatrixXd selectivedynamics;
  if (has_selectivedynamics) {
    AnisoValTraits traits(selectivedynamics_key);
    selectivedynamics =
        mapped_structure.atom_info.properties.at(selectivedynamics_key);
    default_selectivedynamics_atom_value = Eigen::VectorXd::Zero(traits.dim());
  }

  // set occupation - by matching atom_info.names to Molecule::name()
  Index l = 0;
  for (std::string const &name : mapped_structure.atom_info.names) {
    xtal::UnitCellCoord bijk = converter(l);
    xtal::Site const &site = m_xtal_prim->basis()[bijk.sublattice()];
    Index s = 0;
    for (auto const &mol : site.occupant_dof()) {
      if (mol.name() != name) {
        ++s;
        continue;
      }

      // check properties
      auto const &properties = mol.atom(0).properties();

      // selectivedyanmics must match exactly if input
      if (has_selectivedynamics) {
        Eigen::VectorXd ideal_atom_value = default_selectivedynamics_atom_value;
        auto property_it = properties.find(selectivedynamics_key);
        if (property_it != properties.end()) {
          ideal_atom_value = property_it->second.value();
        }

        if (!almost_equal(ideal_atom_value, selectivedynamics.col(l))) {
          ++s;
          continue;
        }
      }

      // name (and selectivedynamics if present ) matches -> break loop
      break;
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
  return default_make_local_dof_values(mapped_structure, m_current_supercell);
}

std::map<std::string, Eigen::MatrixXd>
FromIsotropicAtomicStructure::make_local_properties(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_local_properties(mapped_structure, m_current_supercell);
}

std::map<std::string, Eigen::VectorXd>
FromIsotropicAtomicStructure::make_global_dof_values(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_global_dof_values(mapped_structure);
}

std::map<std::string, Eigen::VectorXd>
FromIsotropicAtomicStructure::make_global_properties(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_global_properties(mapped_structure);
}

/// \brief Constructor
///
/// Requires the following, or throws:
/// - Prim is atomic, with discrete occupants with magspin properties
/// - The mapped_structure has magspin atom properties, of the same flavor
///
/// \param _prim The prim
/// \param _supercells Shared pointer to a SupercellSet. This is an
///     optional method to avoid unnecessarily creating duplicate
///     Supercell. An empty SupercellSet is constructed by default.
/// \param _magspin_tol Maximum allowed difference when mapping magspin,
///     using `(calculated_magspin - ideal_magspin).norm()`.
FromDiscreteMagneticAtomicStructure::FromDiscreteMagneticAtomicStructure(
    std::shared_ptr<Prim const> const &_prim,
    std::shared_ptr<SupercellSet> _supercells, double _magspin_tol)
    : FromStructure(throw_if_equal_to_nullptr(
          _prim,
          "Error in FromDiscreteMagneticAtomicStructure constructor: _prim == "
          "nullptr")),
      m_supercells(_supercells),
      m_magspin_tol(_magspin_tol) {
  if (!m_prim->is_atomic || m_prim->magspin_info.has_continuous_magspin_dof ||
      !m_prim->magspin_info.has_discrete_atomic_magspin_occupants ||
      !m_prim->magspin_info.discrete_atomic_magspin_key.has_value()) {
    throw std::runtime_error(
        "Error in FromDiscreteMagneticAtomicStructure constructor: not a "
        "discrete magnetic atomic prim");
  }

  if (m_supercells == nullptr) {
    m_supercells = std::make_shared<SupercellSet>(m_prim);
  }
}

/// \brief Construct a configuration with properties from a mapped structure
ConfigurationWithProperties FromDiscreteMagneticAtomicStructure::operator()(
    xtal::SimpleStructure const &mapped_structure) {
  // this must be first so other methods can use m_current_supercell
  m_current_supercell =
      m_supercells->insert(this->make_supercell(mapped_structure))
          .first->supercell;

  // Copy <flavor>magspin from atom_info.properties
  std::string magspin_key =
      m_prim->magspin_info.discrete_atomic_magspin_key.value();
  auto magspin_it = mapped_structure.atom_info.properties.find(magspin_key);
  if (magspin_it == mapped_structure.atom_info.properties.end()) {
    std::stringstream msg;
    msg << "Error in FromDiscreteMagneticAtomicStructure: "
        << "mapped_structure does not have atom properties '" << magspin_key
        << "'";
    throw std::runtime_error(msg.str());
  }
  Eigen::MatrixXd magspin = magspin_it->second;

  // Exclude <flavor>magspin from local properties
  std::set<std::string> excluded;
  excluded.insert(magspin_key);

  clexulator::ConfigDoFValues dof_values;
  dof_values.occupation = this->make_occupation(mapped_structure, magspin);
  dof_values.local_dof_values = this->make_local_dof_values(mapped_structure);
  dof_values.global_dof_values = this->make_global_dof_values(mapped_structure);

  return ConfigurationWithProperties(
      Configuration(m_current_supercell, dof_values),
      this->make_local_properties(mapped_structure, excluded),
      this->make_global_properties(mapped_structure));
}

/// \brief Return shared pointer to all supercells
std::shared_ptr<SupercellSet> FromDiscreteMagneticAtomicStructure::supercells()
    const {
  return m_supercells;
}

std::runtime_error FromDiscreteMagneticAtomicStructure::error(
    std::string what) const {
  std::string prefix = "Error in FromIsotropicAtomicStructure: ";
  return std::runtime_error(prefix + what);
}

Eigen::VectorXi FromDiscreteMagneticAtomicStructure::make_occupation(
    xtal::SimpleStructure const &mapped_structure,
    Eigen::MatrixXd const &magspin) const {
  auto const &converter = m_current_supercell->unitcellcoord_index_converter;
  Index n_sites = converter.total_sites();
  validate_atom_names_or_throw(mapped_structure.atom_info.names, n_sites);

  Eigen::VectorXi occupation(n_sites);

  // collect magspin info
  std::string magspin_key =
      m_prim->magspin_info.discrete_atomic_magspin_key.value();
  Eigen::VectorXd default_magspin_atom_value =
      Eigen::VectorXd::Zero(magspin.rows());

  // collect selectivedynamics info
  std::string selectivedynamics_key = "selectivedynamics";
  bool has_selectivedynamics =
      mapped_structure.atom_info.properties.count(selectivedynamics_key);
  Eigen::VectorXd default_selectivedynamics_atom_value;
  Eigen::MatrixXd selectivedynamics;
  if (has_selectivedynamics) {
    AnisoValTraits traits(selectivedynamics_key);
    selectivedynamics =
        mapped_structure.atom_info.properties.at(selectivedynamics_key);
    default_selectivedynamics_atom_value = Eigen::VectorXd::Zero(traits.dim());
  }

  // set occupation - by matching:
  // - atom_info.names to mol.atom(0).name()
  // - if has_selectivedynamics:
  //   - atom_info.properties.at("selectivedynamics") must match ideal value
  // - magspin.col(i) to mol.atom(0).properties().at(magspin_key)
  //   - magnitude of the difference must be less than magspin_tol
  //   - choose closest with matching name
  Index l = 0;
  for (std::string const &name : mapped_structure.atom_info.names) {
    xtal::UnitCellCoord bijk = converter(l);
    xtal::Site const &site = m_xtal_prim->basis()[bijk.sublattice()];
    Index s = 0;  // current occupation index

    // best occupation index (default = not found)
    Index best_s = site.occupant_dof().size();
    double best_magspin_diff = std::numeric_limits<double>::infinity();
    for (auto const &mol : site.occupant_dof()) {
      if (mol.atom(0).name() != name) {
        ++s;
        continue;
      }

      // check properties
      auto const &properties = mol.atom(0).properties();

      // selectivedyanmics must match exactly if input
      if (has_selectivedynamics) {
        Eigen::VectorXd ideal_atom_value = default_selectivedynamics_atom_value;
        auto property_it = properties.find(selectivedynamics_key);
        if (property_it != properties.end()) {
          ideal_atom_value = property_it->second.value();
        }
        if (!almost_equal(ideal_atom_value, selectivedynamics.col(l))) {
          ++s;
          continue;
        }
      }

      // magspin must be within tol, then closest is best
      auto magspin_it = properties.find(magspin_key);
      Eigen::VectorXd ideal_magspin = default_magspin_atom_value;
      if (magspin_it != properties.end()) {
        ideal_magspin = magspin_it->second.value();
      }
      double magspin_diff = (magspin.col(l) - ideal_magspin).norm();

      // do not consider if magspin_diff >= magspin_tol
      if (magspin_diff >= m_magspin_tol) {
        ++s;
        continue;
      }

      // if better diff, store best occupation index
      if (magspin_diff < best_magspin_diff) {
        best_magspin_diff = magspin_diff;
        best_s = s;
      }
      ++s;
    }
    if (best_s == site.occupant_dof().size()) {
      std::stringstream msg;
      msg << "Failed constructing occupation: atom type '" << name
          << "' with magspin == " << magspin.col(l).transpose()
          << " is not allowed on site " << l << ".";
      msg << "Allowed: \n";
      for (auto const &mol : site.occupant_dof()) {
        auto const &properties = mol.atom(0).properties();
        auto magspin_it = properties.find(magspin_key);
        Eigen::VectorXd ideal_magspin = default_magspin_atom_value;
        if (magspin_it != properties.end()) {
          ideal_magspin = magspin_it->second.value();
        }
        msg << "- atom type: " << mol.atom(0).name()
            << " magspin: " << ideal_magspin.transpose() << "\n";
      }
      msg << "\n";
      throw this->error(msg.str());
    }
    occupation(l) = best_s;
    ++l;
  }
  return occupation;
}

std::map<std::string, Eigen::MatrixXd>
FromDiscreteMagneticAtomicStructure::make_local_dof_values(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_local_dof_values(mapped_structure, m_current_supercell);
}

std::map<std::string, Eigen::MatrixXd>
FromDiscreteMagneticAtomicStructure::make_local_properties(
    xtal::SimpleStructure const &mapped_structure,
    std::set<std::string> excluded) const {
  return default_make_local_properties(mapped_structure, m_current_supercell,
                                       excluded);
}

std::map<std::string, Eigen::VectorXd>
FromDiscreteMagneticAtomicStructure::make_global_dof_values(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_global_dof_values(mapped_structure);
}

std::map<std::string, Eigen::VectorXd>
FromDiscreteMagneticAtomicStructure::make_global_properties(
    xtal::SimpleStructure const &mapped_structure) const {
  return default_make_global_properties(mapped_structure);
}

}  // namespace config
}  // namespace CASM
