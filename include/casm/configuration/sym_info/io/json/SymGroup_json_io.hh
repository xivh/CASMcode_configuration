#ifndef CASM_sym_info_SymGroup_json_io
#define CASM_sym_info_SymGroup_json_io

#include <memory>

#include "casm/configuration/sym_info/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;
template <typename T>
struct jsonConstructor;
class jsonParser;

/// \brief Write ClusterSpecs to JSON object
jsonParser &to_json(std::shared_ptr<sym_info::SymGroup const> const &sym_group,
                    jsonParser &json, xtal::Lattice const &lattice);

/// \brief Read from JSON
void from_json(std::shared_ptr<sym_info::SymGroup const> &sym_group,
               jsonParser const &json, xtal::Lattice const &lattice);

template <>
struct jsonConstructor<std::shared_ptr<sym_info::SymGroup const>> {
  /// \brief Construct from JSON
  static std::shared_ptr<sym_info::SymGroup const> from_json(
      jsonParser const &json, xtal::Lattice const &lattice);
};

/// \brief Parse ClusterSpecs from JSON
void parse(InputParser<std::shared_ptr<sym_info::SymGroup const>> &parser,
           xtal::Lattice const &lattice);

}  // namespace CASM

#endif
