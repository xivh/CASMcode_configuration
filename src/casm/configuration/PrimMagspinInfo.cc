#include "casm/configuration/PrimMagspinInfo.hh"

#include "casm/crystallography/BasicStructure.hh"

namespace CASM {
namespace config {

PrimMagspinInfo::PrimMagspinInfo(BasicStructure const &prim) {
  // find continuous magspin flavors
  std::set<std::string> continuous_flavors;
  for (auto const &site : prim.basis()) {
    for (auto const &pair : site.dofs()) {
      std::string key = AnisoValTraits::name_suffix(pair.first);
      if (key.size() != pair.first.size()) {
        std::stringstream msg;
        msg << "Error in PrimMagspinInfo: In a prim, DoF may not be specified "
               "with a qualifier (found '"
            << pair.first << "').";
        throw std::runtime_error(msg.str());
      }
      auto pos = key.find("magspin");
      if (pos != std::string::npos) {
        continuous_flavors.insert(key.substr(0, pos));
      }
    }
  }
  if (continuous_flavors.size() > 1) {
    throw std::runtime_error(
        "Error in PrimMagspinInfo: multiple continuous magspin "
        "flavors found, which is not allowed.");
  } else if (continuous_flavors.size() == 1) {
    std::string flavor = *continuous_flavors.begin();
    has_continuous_magspin_dof = true;
    continuous_magspin_flavor = flavor;
    continuous_magspin_key = flavor + "magspin";
  } else {
    has_continuous_magspin_dof = false;
  }

  // find discrete magspin flavors
  std::vector<xtal::Molecule> molecules = struc_molecule(prim);
  if (molecules.empty()) {
    has_discrete_atomic_magspin_occupants = false;
  } else {
    std::set<std::string> discrete_flavors;
    for (auto const &mol : molecules) {
      if (mol.size() != 1 || !mol.properties().empty()) {
        has_discrete_atomic_magspin_occupants = false;
        break;
      }
      std::map<std::string, xtal::SpeciesProperty> const
          &discrete_atomic_properties = mol.atom(0).properties();
      for (auto const &pair : discrete_atomic_properties) {
        std::string key = AnisoValTraits::name_suffix(pair.first);
        if (key.size() != pair.first.size()) {
          std::stringstream msg;
          msg << "Error in PrimMagspinInfo: In a prim, occupants may not be "
                 "defined "
                 "using properties with a qualifier (found '"
              << pair.first << "').";
          throw std::runtime_error(msg.str());
        }
        auto pos = key.find("magspin");
        if (pos != std::string::npos) {
          discrete_flavors.insert(key.substr(0, pos));
        }
      }
    }
    if (discrete_flavors.size() > 1) {
      throw std::runtime_error(
          "Error in PrimMagspinInfo: multiple discrete magspin "
          "flavors found, which is not allowed.");
    } else if (discrete_flavors.size() == 1) {
      std::string flavor = *discrete_flavors.begin();
      has_discrete_atomic_magspin_occupants = true;
      discrete_atomic_magspin_flavor = flavor;
      discrete_atomic_magspin_key = flavor + "magspin";
    } else {
      has_discrete_atomic_magspin_occupants = false;
    }
  }
}

}  // namespace config
}  // namespace CASM
