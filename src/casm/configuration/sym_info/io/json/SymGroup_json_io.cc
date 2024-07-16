#include "casm/configuration/sym_info/io/json/SymGroup_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/sym_info/factor_group.hh"
#include "casm/crystallography/SymInfo.hh"
#include "casm/crystallography/io/SymInfo_json_io.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

namespace {

void _write_symop(sym_info::SymGroup const &grp, Index i, jsonParser &j,
                  xtal::Lattice const &lattice,
                  std::vector<Index> const &conjugacy_class_of_op,
                  std::vector<xtal::SymInfo> const &sym_info_of_op) {
  j = jsonParser::object();

  const xtal::SymOp &op = grp.element[i];

  j["master_group_index"] =
      grp.head_group_index[i];  // op.master_group_index();

  to_json(op.matrix, j["CART"]["matrix"]);
  to_json_array(op.translation, j["CART"]["tau"]);
  to_json(op.is_time_reversal_active, j["CART"]["time_reversal"]);

  to_json(lattice.inv_lat_column_mat() * op.matrix * lattice.lat_column_mat(),
          j["FRAC"]["matrix"]);
  to_json_array(lattice.inv_lat_column_mat() * op.translation,
                j["FRAC"]["tau"]);
  to_json(op.is_time_reversal_active, j["FRAC"]["time_reversal"]);

  to_json(conjugacy_class_of_op[i] + 1, j["info"]["conjugacy_class"]);
  to_json(grp.inv(i) + 1, j["info"]["inverse_operation"]);

  to_json(sym_info_of_op[i], j["info"]);

  to_json(to_brief_unicode(sym_info_of_op[i], xtal::SymInfoOptions(CART)),
          j["info"]["brief"]["CART"]);
  to_json(to_brief_unicode(sym_info_of_op[i], xtal::SymInfoOptions(FRAC)),
          j["info"]["brief"]["FRAC"]);
}

}  // namespace

/// \brief Write ClusterSpecs to JSON object
jsonParser &to_json(std::shared_ptr<sym_info::SymGroup const> const &sym_group,
                    jsonParser &json, xtal::Lattice const &lattice) {
  std::vector<std::vector<Index>> conjugacy_classes =
      make_conjugacy_classes(*sym_group);

  std::vector<Index> conjugacy_class_of_op(sym_group->element.size(), -1);
  for (Index c = 0; c < conjugacy_classes.size(); ++c) {
    for (Index op_index : conjugacy_classes[c]) {
      conjugacy_class_of_op[op_index] = c;
    }
  }

  std::vector<xtal::SymInfo> sym_info_of_op;
  sym_info_of_op.reserve(sym_group->element.size());
  for (auto const &op : sym_group->element) {
    sym_info_of_op.emplace_back(op, lattice);
  }

  json = jsonParser::object();

  {
    jsonParser &json_ops = json["group_operations"];
    for (int i = 0; i < sym_group->element.size(); i++) {
      std::string op_name =
          "op_" + to_sequential_string(i + 1, sym_group->element.size());
      _write_symop(*sym_group, i, json_ops[op_name], lattice,
                   conjugacy_class_of_op, sym_info_of_op);
    }
  }

  {
    jsonParser &json_info = json["group_classification"];
    // json_info["name"].put_null();
    // json_info["latex_name"].put_null();
    // json_info["periodicity"] = "PERIODIC";
  }

  {
    jsonParser &json_struc = json["group_structure"];

    for (Index c = 0; c < conjugacy_classes.size(); ++c) {
      std::string class_name =
          "class_" + to_sequential_string(c + 1, conjugacy_classes.size());
      jsonParser &json_class = json_struc["conjugacy_classes"][class_name];

      xtal::SymInfo const &info = sym_info_of_op[conjugacy_classes[c][0]];
      json_class["operation_type"] = info.op_type;
      if (info.op_type == xtal::symmetry_type::rotation_op ||
          info.op_type == xtal::symmetry_type::screw_op ||
          info.op_type == xtal::symmetry_type::rotoinversion_op) {
        json_class["rotation_angle"] = info.angle;
      }
      json_class["operations"].put_array();
      for (Index o : conjugacy_classes[c]) {
        json_class["operations"].push_back(o + 1);
      }
    }

    json_struc["multiplication_table"] = sym_group->multiplication_table;
  }

  return json;
}

/// \brief Read from JSON
void from_json(std::shared_ptr<sym_info::SymGroup const> &sym_group,
               jsonParser const &json, xtal::Lattice const &lattice) {
  sym_group =
      jsonConstructor<std::shared_ptr<sym_info::SymGroup const>>::from_json(
          json, lattice);
}

/// \brief Construct from JSON
std::shared_ptr<sym_info::SymGroup const>
jsonConstructor<std::shared_ptr<sym_info::SymGroup const>>::from_json(
    jsonParser const &json, xtal::Lattice const &lattice) {
  InputParser<std::shared_ptr<sym_info::SymGroup const>> parser{json, lattice};
  std::stringstream ss;
  ss << "Error: Invalid SymGroup JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Parse ClusterSpecs from JSON
void parse(InputParser<std::shared_ptr<sym_info::SymGroup const>> &parser,
           xtal::Lattice const &lattice) {
  if (!parser.self.contains("group_operations")) {
    std::stringstream ss;
    ss << "Error reading SymGroup from JSON: missing group_operations";
    parser.insert_error("group_operations", ss.str());
  }

  std::vector<xtal::SymOp> elements;
  Index n_elements = parser.self["group_operations"].size();
  for (Index i = 0; i < parser.self["group_operations"].size(); ++i) {
    fs::path op_path{"group_operations"};
    std::string op_name = "op_" + to_sequential_string(i + 1, n_elements);

    Eigen::Matrix3d matrix;
    parser.require(matrix, op_path / op_name / "CART" / "matrix");

    Eigen::Vector3d translation;
    parser.require(translation, op_path / op_name / "CART" / "tau");

    bool is_time_reversal_active;
    parser.require(is_time_reversal_active,
                   op_path / op_name / "CART" / "time_reversal");

    elements.emplace_back(matrix, translation, is_time_reversal_active);
  }

  if (parser.valid()) {
    parser.value = std::make_unique<std::shared_ptr<sym_info::SymGroup const>>(
        sym_info::make_symgroup(elements, lattice));
  }
}

}  // namespace CASM
