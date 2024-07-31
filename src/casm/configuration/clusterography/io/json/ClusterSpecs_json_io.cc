#include "casm/configuration/clusterography/io/json/ClusterSpecs_json_io.hh"

#include <optional>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/configuration/clusterography/ClusterInvariants.hh"
#include "casm/configuration/clusterography/ClusterSpecs.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/io/json/IntegralClusterOrbitGenerator_json_io.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/clusterography/orbits.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/global/enum/json_io.hh"

namespace CASM {

/// \brief Write ClusterSpecs to JSON object
jsonParser &to_json(clust::ClusterSpecs const &cluster_specs,
                    jsonParser &json) {
  // generating_group
  json["generating_group"] = cluster_specs.generating_group->head_group_index;

  // site_filter:
  json["site_filter_method"] = cluster_specs.site_filter_method;

  // orbit_branch_specs
  json["orbit_branch_specs"];
  for (Index i = 0; i < cluster_specs.max_length.size(); ++i) {
    std::string branch = std::to_string(i);
    json["orbit_branch_specs"][branch]["max_length"] =
        cluster_specs.max_length[i];
    if (i < cluster_specs.cutoff_radius.size()) {
      json["orbit_branch_specs"][branch]["cutoff_radius"] =
          cluster_specs.cutoff_radius[i];
    }
  }

  auto const &prim = *cluster_specs.prim;

  // orbit_specs
  if (cluster_specs.custom_generators.size()) {
    to_json(cluster_specs.custom_generators, json["orbit_specs"], prim);
  }

  // phenomenal
  if (cluster_specs.phenomenal.has_value()) {
    to_json(*cluster_specs.phenomenal, json["phenomenal"], prim);
  }

  // include_phenomenal_sites
  if (cluster_specs.include_phenomenal_sites) {
    to_json(cluster_specs.include_phenomenal_sites,
            json["include_phenomenal_sites"]);
  }
  return json;
}

/// \brief Read from JSON
void from_json(
    clust::ClusterSpecs &cluster_specs, jsonParser const &json,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep) {
  cluster_specs = jsonConstructor<clust::ClusterSpecs>::from_json(
      json, prim, prim_factor_group, unitcellcoord_symgroup_rep);
}

/// \brief Construct from JSON
clust::ClusterSpecs jsonConstructor<clust::ClusterSpecs>::from_json(
    jsonParser const &json,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep) {
  InputParser<clust::ClusterSpecs> parser{json, prim, prim_factor_group,
                                          unitcellcoord_symgroup_rep};
  std::stringstream ss;
  ss << "Error: Invalid cluster_specs JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

/// \brief Parse ClusterSpecs from JSON
void parse(
    InputParser<clust::ClusterSpecs> &parser,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::shared_ptr<clust::SymGroup const> const &prim_factor_group,
    std::vector<xtal::UnitCellCoordRep> const &unitcellcoord_symgroup_rep) {
  // get phenomenal (optional)
  std::optional<clust::IntegralCluster> phenomenal;
  parser.optional(phenomenal, "phenomenal", *prim);

  // get include_phenomenal_sites (optional)
  bool include_phenomenal_sites = false;
  parser.optional(include_phenomenal_sites, "include_phenomenal_sites");

  // get generating_group (optional)
  std::shared_ptr<clust::SymGroup const> generating_group;
  std::set<Index> head_group_index;
  if (parser.self.contains("generating_group")) {
    parser.optional(head_group_index, "generating_group");
  } else {
    for (Index i : prim_factor_group->head_group_index) {
      head_group_index.insert(i);
    }
  }
  if (phenomenal.has_value()) {
    // need to get elements consistent with phenomenal cluster invariance
    std::vector<xtal::SymOp> cluster_group_elements;
    for (Index i : head_group_index) {
      cluster_group_elements.push_back(clust::make_cluster_group_element(
          *phenomenal, prim->lattice().lat_column_mat(),
          prim_factor_group->element[i], unitcellcoord_symgroup_rep[i]));
    }
    generating_group = std::make_shared<clust::SymGroup const>(
        prim_factor_group, cluster_group_elements, head_group_index);
  } else {
    std::vector<xtal::SymOp> element;
    for (Index i : head_group_index) {
      element.push_back(prim_factor_group->element[i]);
    }
    generating_group = std::make_shared<clust::SymGroup const>(
        prim_factor_group, element, head_group_index);
  }

  // site_filter:
  std::string site_filter_method = "dof_sites";
  parser.optional<std::string>(site_filter_method,
                               fs::path("site_filter_method"));

  clust::SiteFilterFunction site_filter;
  if (site_filter_method == "dof_sites") {
    site_filter = clust::dof_sites_filter();
  } else if (site_filter_method == "alloy_sites") {
    site_filter = clust::alloy_sites_filter;
  } else if (site_filter_method == "all_sites") {
    site_filter = clust::all_sites_filter;
  } else {
    std::stringstream ss;
    ss << "Error reading ClusterSpecs from JSON: site_filter_method="
       << site_filter_method << " is not recognized";
    parser.insert_error("site_filter_method", ss.str());
  }

  // orbit_branch_specs
  std::vector<double> max_length;
  std::vector<double> cutoff_radius;
  if (parser.self.contains("orbit_branch_specs")) {
    bool has_cutoff_radius = false;
    jsonParser const &j = parser.self["orbit_branch_specs"];
    std::size_t pos{};
    for (auto it = j.begin(); it != j.end(); ++it) {
      try {
        const int i{std::stoi(it.name(), &pos)};
        while (max_length.size() < i + 1) {
          max_length.push_back(0.0);
          cutoff_radius.push_back(0.0);
        }

        if (!it->is_obj()) {
          continue;
        }
        if (it->contains("cutoff_radius")) {
          has_cutoff_radius = true;
        }

        // read "max_length" & "cutoff_radius" (defaults = 0.0)
        fs::path option = fs::path("orbit_branch_specs") / it.name();
        double _max_length = 0.0;
        parser.optional(_max_length, option / "max_length");
        max_length[i] = _max_length;

        double _cutoff_radius = 0.0;
        parser.optional(_cutoff_radius, option / "cutoff_radius");
        cutoff_radius[i] = _cutoff_radius;
      } catch (std::invalid_argument const &ex) {
        continue;
      }
    }

    // if "cutoff_radius" never included in input, clear
    if (!has_cutoff_radius) {
      cutoff_radius.clear();
    }
  } else if (all_local_dof_types(*prim).size()) {
    std::stringstream ss;
    ss << "Warning reading ClusterSpecs from JSON: "
       << "prim has local DoF and `orbit_branch_specs` does not exist. "
       << "To suppress this warning, include `\"orbit_branch_specs\": {}`";
    parser.insert_warning("orbit_branch_specs", ss.str());
  }

  // orbit_specs (optional)
  // empty by default
  std::vector<clust::IntegralClusterOrbitGenerator> default_custom_generators;
  auto custom_generators_parser =
      parser.subparse_else<std::vector<clust::IntegralClusterOrbitGenerator>>(
          "orbit_specs", default_custom_generators, *prim);

  if (!parser.valid()) {
    return;
  }
  parser.value =
      notstd::make_unique<clust::ClusterSpecs>(prim, generating_group);

  parser.value->site_filter_method = site_filter_method;
  parser.value->site_filter = site_filter;
  parser.value->max_length = max_length;
  parser.value->custom_generators = *custom_generators_parser->value;

  // --- Local clusters only ---
  parser.value->phenomenal = phenomenal;
  parser.value->include_phenomenal_sites = include_phenomenal_sites;
  parser.value->cutoff_radius = cutoff_radius;
}

}  // namespace CASM
