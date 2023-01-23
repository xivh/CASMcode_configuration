#ifndef CASM_config_make_simple_structure
#define CASM_config_make_simple_structure

#include <map>
#include <string>

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
}

namespace config {

struct Configuration;

/// \brief Convert a Configuration to a SimpleStructure
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration,
    std::map<std::string, Eigen::MatrixXd> const &local_properties = {},
    std::map<std::string, Eigen::VectorXd> const &global_properties = {});

}  // namespace config
}  // namespace CASM

#endif
