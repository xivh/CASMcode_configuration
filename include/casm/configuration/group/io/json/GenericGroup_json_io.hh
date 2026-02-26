#ifndef CASM_group_GenericGroup_json_io
#define CASM_group_GenericGroup_json_io

#include <memory>

#include "casm/configuration/group/definitions.hh"

namespace CASM {

template <typename T>
class InputParser;
template <typename T>
struct jsonConstructor;
class jsonParser;

namespace group {

struct GenericGroup;
struct Subset;

}  // namespace group

/// \brief Write GenericGroup to JSON object
jsonParser &to_json(std::shared_ptr<group::GenericGroup const> const &group,
                    jsonParser &json);

/// \brief Read from JSON
void from_json(std::shared_ptr<group::GenericGroup const> &group,
               jsonParser const &json);

template <>
struct jsonConstructor<std::shared_ptr<group::GenericGroup const>> {
  /// \brief Construct from JSON
  static std::shared_ptr<group::GenericGroup const> from_json(
      jsonParser const &json);
};

/// \brief Parse GenericGroup from JSON
void parse(InputParser<std::shared_ptr<group::GenericGroup const>> &parser);

}  // namespace CASM

#endif
