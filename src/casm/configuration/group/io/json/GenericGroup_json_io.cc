#include "casm/configuration/group/io/json/GenericGroup_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/subgroups.hh"

namespace CASM {

/// \brief Write GenericGroup to JSON object
jsonParser &to_json(std::shared_ptr<group::GenericGroup const> const &group,
                    jsonParser &json) {
  std::vector<std::vector<Index>> conjugacy_classes =
      make_conjugacy_classes(*group);

  json = jsonParser::object();

  {
    jsonParser &json_struc = json["group_structure"];

    for (Index c = 0; c < conjugacy_classes.size(); ++c) {
      std::string class_name =
          "class_" + to_sequential_string(c + 1, conjugacy_classes.size());
      jsonParser &json_class = json_struc["conjugacy_classes"][class_name];

      json_class["operations"].put_array();
      for (Index o : conjugacy_classes[c]) {
        json_class["operations"].push_back(o + 1);
      }
    }

    json_struc["multiplication_table"] = group->multiplication_table;
  }

  return json;
}

/// \brief Read from JSON
void from_json(std::shared_ptr<group::GenericGroup const> &group,
               jsonParser const &json) {
  group =
      jsonConstructor<std::shared_ptr<group::GenericGroup const>>::from_json(
          json);
}

/// \brief Construct from JSON
std::shared_ptr<group::GenericGroup const>
jsonConstructor<std::shared_ptr<group::GenericGroup const>>::from_json(
    jsonParser const &json) {
  InputParser<std::shared_ptr<group::GenericGroup const>> parser{json};
  std::stringstream ss;
  ss << "Error: Invalid GenericGroup JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Parse GenericGroup from JSON
void parse(InputParser<std::shared_ptr<group::GenericGroup const>> &parser) {
  fs::path group_structure_path = fs::path{"group_structure"};
  fs::path mtable_path = group_structure_path / "multiplication_table";
  if (parser.self.find(mtable_path) == parser.self.end()) {
    std::stringstream ss;
    ss << "Error reading GenericGroup from JSON: missing multiplication_table";
    parser.insert_error(mtable_path, ss.str());
  }

  std::vector<std::vector<Index>> multiplication_table;
  parser.require(multiplication_table, mtable_path);

  if (parser.valid()) {
    parser.value = std::make_unique<std::shared_ptr<group::GenericGroup const>>(
        std::make_shared<group::GenericGroup>(multiplication_table));
  }
}

}  // namespace CASM
