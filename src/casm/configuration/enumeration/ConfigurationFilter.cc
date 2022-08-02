#include "casm/configuration/enumeration/ConfigurationFilter.hh"

#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"
#include "casm/configuration/copy_configuration.hh"

namespace CASM {
namespace config {

bool UniqueConfigurationFilter::operator()(
    Configuration const &configuration) const {
  if (!is_primitive(configuration)) {
    return false;
  }
  auto begin = SupercellSymOp::begin(configuration.supercell);
  auto end = SupercellSymOp::end(configuration.supercell);
  if (!is_canonical(configuration, begin, end)) {
    return false;
  }
  return true;
}

bool GenericConfigurationFilter::operator()(
    Configuration const &configuration) const {
  if (primitive_only && !is_primitive(configuration)) {
    return false;
  }
  auto begin = SupercellSymOp::begin(configuration.supercell);
  auto end = SupercellSymOp::end(configuration.supercell);
  if (canonical_only && !is_canonical(configuration, begin, end)) {
    return false;
  }
  return f(configuration);
}

}  // namespace config
}  // namespace CASM
