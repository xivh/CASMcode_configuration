#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/clusterography/ClusterInvariants.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/global/enum/json_io.hh"

namespace CASM {

namespace {

clust::ClusterInvariants _make_cluster_invariants(
    clust::IntegralCluster const &cluster, xtal::BasicStructure const &prim,
    std::optional<clust::IntegralCluster> phenomenal) {
  if (phenomenal) {
    return clust::ClusterInvariants(cluster, *phenomenal, prim);
  } else {
    return clust::ClusterInvariants(cluster, prim);
  }
}

}  // anonymous namespace

/// \brief Write IntegralCluster to JSON object
///
/// Format:
/// \code
/// {
///   "min_length" : number,
///   "max_length" : number,
///   "sites" : [
///     [b, i, j, k],
///     ...
///   ]
/// }
/// \endcode
jsonParser &to_json(clust::IntegralCluster const &clust, jsonParser &json,
                    xtal::BasicStructure const &prim,
                    std::optional<clust::IntegralCluster> phenomenal) {
  clust::ClusterInvariants invariants =
      _make_cluster_invariants(clust, prim, phenomenal);
  if (invariants.distances().size()) {
    json["min_length"] = invariants.distances().front();
    json["max_length"] = invariants.distances().back();
  } else {
    json["min_length"] = 0.0;
    json["max_length"] = 0.0;
  }
  json["distances"] = invariants.distances();
  if (invariants.phenomenal_distances().size()) {
    json["phenomenal_distances"] = invariants.phenomenal_distances();
  }
  json["sites"].put_array(clust.begin(), clust.end());
  return json;
}

/// \brief Read IntegralCluster from JSON
///
/// Format:
/// \code
/// {
///   "min_length" : number,
///   "max_length" : number,
///   "coordinate_mode" : ("FRAC", "CART", "INT" (default)) (optional)
///   "sites" : [
///     [b, i, j, k],
///     ...
///   ]
/// }
/// \endcode
///
/// - Also accepts "prototype" in place of "sites"
void from_json(clust::IntegralCluster &clust, const jsonParser &json,
               xtal::BasicStructure const &prim) {
  clust = jsonConstructor<clust::IntegralCluster>::from_json(json, prim);
}

clust::IntegralCluster jsonConstructor<clust::IntegralCluster>::from_json(
    const jsonParser &json, xtal::BasicStructure const &prim) {
  InputParser<clust::IntegralCluster> parser{json, prim};
  std::stringstream ss;
  ss << "Error: Invalid cluster JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Parse IntegralCluster from JSON
///
/// Format:
/// \code
/// {
///   "coordinate_mode" : ("FRAC", "CART", "INT" (default)) (optional)
///   "sites" : [
///     [b, i, j, k],
///     ...
///   ]
/// }
/// \endcode
///
/// - Also accepts "prototype" in place of "sites"
void parse(InputParser<clust::IntegralCluster> &parser,
           xtal::BasicStructure const &prim) {
  std::string name;
  const jsonParser &json = parser.self;
  if (json.contains("sites")) {
    name = "sites";
  } else if (json.contains("prototype")) {
    name = "prototype";
  } else {
    parser.error.insert(
        "Error reading IntegralCluster from JSON: Expected 'sites' or "
        "'prototype' containing a list of cluster site coordinates.");
    return;
  }

  double xtal_tol = prim.lattice().tol();
  CASM::COORD_TYPE coord_type;
  parser.optional_else(coord_type, "coordinate_mode", INTEGRAL);
  if (!parser.valid()) {
    return;
  }

  parser.value = notstd::make_unique<clust::IntegralCluster>();
  auto &clust = *parser.value;

  if (coord_type == INTEGRAL) {
    parser.require(clust.elements(), name);
  } else {
    std::vector<Eigen::VectorXd> coord_vec;
    parser.require(coord_vec, name);

    try {
      for (const auto &coord : coord_vec) {
        if (coord.size() != 3) {
          parser.error.insert(
              "Error: cluster site coordinates have wrong dimension");
          parser.value.reset();
          return;
        }
        xtal::Coordinate tcoord{coord, prim.lattice(), coord_type};
        clust.elements().emplace_back(
            xtal::UnitCellCoord::from_coordinate(prim, tcoord, xtal_tol));
      }
    } catch (std::exception &e) {
      parser.error.insert("Error: could not read coordinates from '" + name +
                          "'");
      return;
    }
  }
  return;
}
}  // namespace CASM
