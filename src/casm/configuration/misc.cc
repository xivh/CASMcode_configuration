#include "casm/configuration/misc.hh"

#include "casm/casm_io/Log.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
namespace config {

AnisoValTraits get_local_traits_or_throw(std::string key) {
  try {
    return AnisoValTraits(key);
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "CASM does not know how to transform the local property '" << key
        << "'. The property name suffix must be the name of a local property "
           "that CASM can transform.";
    CASM::err_log() << std::endl;
    CASM::err_log().paragraph(msg.str());
    CASM::err_log() << std::endl;
    CASM::err_log() << "Local properties include: " << std::endl;
    for (auto const &pair : AnisoValTraits::registered()) {
      std::string name = pair.first;
      auto const &traits = pair.second;
      if (traits.global()) continue;
      CASM::err_log() << "- " << name << std::endl;
    }
    CASM::err_log() << std::endl;
    throw std::runtime_error(std::string("Cannot transform local property '") +
                             key + "'");
  }
}

AnisoValTraits get_global_traits_or_throw(std::string key) {
  try {
    return AnisoValTraits(key);
  } catch (std::exception &e) {
    std::stringstream msg;
    msg << "CASM does not know how to transform the global property '" << key
        << "'. The property name suffix must be the name of a global property "
           "that CASM can transform.";
    CASM::err_log() << std::endl;
    CASM::err_log().paragraph(msg.str());
    CASM::err_log() << std::endl;
    CASM::err_log() << "Global properties include: " << std::endl;
    for (auto const &pair : AnisoValTraits::registered()) {
      std::string name = pair.first;
      auto const &traits = pair.second;
      if (!traits.global()) continue;
      CASM::err_log() << "- " << name << std::endl;
    }
    CASM::err_log() << std::endl;
    throw std::runtime_error(std::string("Cannot transform global property '") +
                             key + "'");
  }
}

}  // namespace config
}  // namespace CASM
