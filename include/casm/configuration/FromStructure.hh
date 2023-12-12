#ifndef CASM_config_FromStructure
#define CASM_config_FromStructure

#include <map>
#include <set>
#include <string>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/global/filesystem.hh"

namespace CASM {
class AnisoValTraits;

namespace xtal {
class BasicStructure;
class SimpleStructure;
}  // namespace xtal

namespace config {
struct ConfigurationWithProperties;
struct Prim;
struct Supercell;
struct SupercellRecord;
class SupercellSet;

/// \brief Base class to help construct a configuration from
///     a mapped structure
///
/// Notes:
/// - The `mapped_structure` this class acts on is an structure
///   with only isotropic atoms on the basis sites which has been
///   "mapped", meaning that atoms and atom properties are in
///   the same order as the supercell sites they are mapped to and
///   all properties (including the coordinates, displacements,
///   and lattice vectors) are rotated so that they can be directly
///   copied to configuration DoF or properties.
/// - If mapped_structure has a global strain property (i.e.
///   "Ustrain", "Hstrain", "GLstrain", etc.) it must be of a type
///   recognized by CASM; it is interpreted as the strain applied to
///   the ideal lattice vectors to result in the structure's lattice
///   vectors; and it is converted to "Ustrain".
/// - Solves L_mapped = Ustrain * L_ideal for L_ideal, which must
///   satisfy L_ideal == L_prim * T, where T is an integer matrix.
///
class FromStructure {
 public:
  FromStructure(std::shared_ptr<Prim const> const &_prim);

 protected:
  virtual std::runtime_error error(std::string what) const;

  AnisoValTraits get_local_traits_or_throw(std::string key) const;

  AnisoValTraits get_global_traits_or_throw(std::string key) const;

  Eigen::VectorXd get_Ustrain_vector(
      xtal::SimpleStructure const &mapped_structure) const;

  Eigen::Matrix3d get_Ustrain_matrix(
      xtal::SimpleStructure const &mapped_structure) const;

  std::shared_ptr<Supercell const> make_supercell(
      xtal::SimpleStructure const &mapped_structure) const;

  Eigen::MatrixXd const &get_local_property_or_throw(
      std::string key,
      std::map<std::string, Eigen::MatrixXd> const &properties) const;

  Eigen::MatrixXd const &get_global_property_or_throw(
      std::string key,
      std::map<std::string, Eigen::MatrixXd> const &properties) const;

  void validate_atom_coords_or_throw(
      xtal::SimpleStructure const &mapped_structure,
      std::shared_ptr<Supercell const> const &supercell) const;

  void validate_atom_names_or_throw(std::vector<std::string> const &names,
                                    Index n_sites) const;

  void validate_local_property_or_throw(
      std::string key, Eigen::MatrixXd const &standard_dof_values,
      Index n_sites) const;

  void validate_global_property_or_throw(
      std::string key, Eigen::MatrixXd const &standard_dof_values) const;

  fs::path error_filename() const;

 protected:
  std::shared_ptr<Prim const> m_prim;
  std::shared_ptr<xtal::BasicStructure const> m_xtal_prim;
};

/// \brief Construct a configuration with properties from
///     a mapped structure of isotropic atoms
///
/// Notes:
/// - The `mapped_structure` this class acts on is an structure
///   with only isotropic atoms on the basis sites which has been
///   "mapped", meaning that atoms and atom properties are in
///   the same order as the supercell sites they are mapped to and
///   all properties (including the coordinates, displacements,
///   and lattice vectors) are rotated so that they can be directly
///   copied to configuration DoF or properties.
/// - If mapped_structure has a global strain property (i.e.
///   "Ustrain", "Hstrain", "GLstrain", etc.) it must be of a type
///   recognized by CASM; it is interpreted as the strain applied to
///   the ideal lattice vectors to result in the structure's lattice
///   vectors; and it is converted to "Ustrain".
/// - Solves L_mapped = Ustrain * L_ideal for L_ideal, which must
///   satisfy L_ideal == L_prim * T, where T is an integer matrix.
/// - Atomic "disp" properties are the displacements from the ideal
///   site coordinates to the structure's atomic coordinates.
/// - If a strain DoF exists, then the strain property of the mapped
///   structure is converted to the DoF strain metric.
/// - If a strain property does not exist, then it is assumed that
///   there is no strain and the mapped_structure's lattice is ideal.
///
class FromIsotropicAtomicStructure : public FromStructure {
 public:
  /// \brief Constructor
  FromIsotropicAtomicStructure(
      std::shared_ptr<Prim const> const &_prim,
      std::shared_ptr<SupercellSet> _supercells = nullptr);

  /// \brief Construct a configuration with properties from a mapped structure
  ConfigurationWithProperties operator()(
      xtal::SimpleStructure const &mapped_structure);

  /// \brief Return shared pointer to all supercells
  std::shared_ptr<SupercellSet> supercells() const;

 protected:
  std::runtime_error error(std::string what) const override;

  Eigen::VectorXi make_occupation(
      xtal::SimpleStructure const &mapped_structure) const;

  std::map<std::string, Eigen::MatrixXd> make_local_dof_values(
      xtal::SimpleStructure const &mapped_structure) const;

  std::map<std::string, Eigen::MatrixXd> make_local_properties(
      xtal::SimpleStructure const &mapped_structure) const;

  std::map<std::string, Eigen::VectorXd> make_global_dof_values(
      xtal::SimpleStructure const &mapped_structure) const;

  std::map<std::string, Eigen::VectorXd> make_global_properties(
      xtal::SimpleStructure const &mapped_structure) const;

 private:
  // Current supercell
  std::shared_ptr<Supercell const> m_current_supercell;

  // Set of supercells to avoid unnecessary duplication
  std::shared_ptr<SupercellSet> m_supercells;
};

}  // namespace config
}  // namespace CASM

#endif
