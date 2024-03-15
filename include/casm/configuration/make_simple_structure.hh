#ifndef CASM_config_make_simple_structure
#define CASM_config_make_simple_structure

#include <map>
#include <set>
#include <string>

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
}

namespace config {

struct Configuration;
struct ConfigurationWithProperties;

/// \brief (deprecated) Convert a Configuration to a SimpleStructure
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration,
    std::map<std::string, Eigen::MatrixXd> const &local_properties = {},
    std::map<std::string, Eigen::VectorXd> const &global_properties = {},
    std::string atom_type_naming_method = "chemical_name",
    std::set<std::string> excluded_species = {"Va", "VA", "va"});

/// \brief Construct a SimpleStructure from a configuration of a Prim with
///     atomic occupants
///
/// Notes:
/// - This is a simpler method which only accepts configurations from
///   prim with only atomic occupants
/// - If "disp" is a DoF, it is included in SimpleStructure::atom_info.coords,
///   but not included in SimpleStructure::atom_info.properties.
/// - If "disp" is included in `local_properties`, it is included in both
///   SimpleStructure::atom_info.coords and
///   SimpleStructure::atom_info.properties.
/// - If strain is a DoF, it is included in SimpleStructure::atom_info.coords
///   and SimpleStructure::lat_column_mat, but not included in
///   SimpleStructure::properties
/// - If "Ustrain" is included in `global_properties`, it is included in
///   SimpleStructure::lat_column_mat.
/// - Other DoF are copied to SimpleStructure::atom_info.properties or
///   SimpleStructure::properties.
/// - Deformation gradient applied to ideal lattice vectors is obtained from:
///   1) a strain DoF (this is not copied to structure.properties)
///   2) "Ustrain" in global_properties (this is copied to structure.properties)
///   3) the identify matrix, I
/// - Options for `atom_type_naming_method` are:
///   - "chemical_name" (default): use ``mol.atom(0).name()``, where `mol`
///     is the `xtal::Molecule` occupying the site
///   - "orientation_name": use ``unique_names[b][s]``, where `unique_names` is
///     the unique occupant names obtained from
///     ``xtal::BasicStructure::unique_names``, `b` is the sublattice index of
///     the site, and `s` is the occupation index.
/// - Magnetic spin may be included as a continuous DoF or as discrete DoF
///   by being an property of atomic occupant, but not both.
/// - All magnetic spins must be of the same flavor.
/// - Magnetic spins are included as `local_properties`. Any occupants
//    that do not have a magnetic spin atomic property are
//    given a value of `[0]` (for collinear flavors) or `[0, 0, 0]`
//    (for non-collinear flavors). and set to [0.]
class ToAtomicStructure {
 public:
  /// \brief Constructor
  ToAtomicStructure(std::string atom_type_naming_method = "chemical_name",
                    std::set<std::string> excluded_species = {"Va", "VA",
                                                              "va"});

  /// \brief Convert a ConfigurationWithProperties to a SimpleStructure
  xtal::SimpleStructure operator()(
      ConfigurationWithProperties const &configuration_with_properties);

  /// \brief Convert a Configuration to a SimpleStructure
  xtal::SimpleStructure operator()(
      Configuration const &configuration,
      std::map<std::string, Eigen::MatrixXd> const &local_properties = {},
      std::map<std::string, Eigen::VectorXd> const &global_properties = {});

 private:
  std::string m_atom_type_naming_method;
  std::set<std::string> m_excluded_species;
};

}  // namespace config
}  // namespace CASM

#endif
