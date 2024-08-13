#include "casm/configuration/occ_events/io/json/OccEvent_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/configuration/group/Group.hh"
#include "casm/configuration/group/orbits.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccEventInvariants.hh"
#include "casm/configuration/occ_events/OccEventRep.hh"
#include "casm/configuration/occ_events/OccPosition.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/configuration/occ_events/OccTrajectory.hh"
#include "casm/configuration/occ_events/orbits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {

namespace {  // anonymous for implementation methods

void parse_reservoir_position(InputParser<occ_events::OccPosition> &parser,
                              occ_events::OccSystem const &system) {
  // parse "occupant_index"
  Index occupant_index;
  parser.require(occupant_index, "occupant_index");
  if (!parser.valid()) {
    return;
  }
  std::unique_ptr<occ_events::OccPosition> p;
  try {
    p = std::make_unique<occ_events::OccPosition>(
        system.make_molecule_in_reservoir_position(occupant_index));
  } catch (std::exception &e) {
    parser.insert_error("occupant_index", e.what());
    return;
  }
  if (parser.valid()) {
    parser.value = std::move(p);
  }
}

void parse_non_reservoir_position(InputParser<occ_events::OccPosition> &parser,
                                  occ_events::OccSystem const &system) {
  // parse and validate "coordinate"
  xtal::UnitCellCoord integral_site_coordinate{};
  parser.require(integral_site_coordinate, "coordinate");
  Index b = integral_site_coordinate.sublattice();
  if (b < 0 || b >= system.prim->basis().size()) {
    parser.error.insert("Error: Invalid coordinate");
    return;
  }

  // parse "occupant_index"
  Index occupant_index;
  parser.require(occupant_index, "occupant_index");
  if (occupant_index < 0 ||
      occupant_index >= system.atom_position_to_name_index[b].size()) {
    parser.error.insert("Error: Invalid occupant_index");
    return;
  }
  Index mol_size = system.atom_position_to_name_index[b][occupant_index].size();

  std::unique_ptr<occ_events::OccPosition> p;
  try {
    if (parser.self.contains("molecule")) {
      p = std::make_unique<occ_events::OccPosition>(
          system.make_molecule_position(integral_site_coordinate,
                                        occupant_index));
    } else if (mol_size == 1) {
      p = std::make_unique<occ_events::OccPosition>(
          system.make_atom_position(integral_site_coordinate, occupant_index));
    } else {
      // parse "atom_position_index"
      Index atom_position_index = -1;
      parser.require(atom_position_index, "atom_position_index");
      p = std::make_unique<occ_events::OccPosition>(system.make_atom_position(
          integral_site_coordinate, occupant_index, atom_position_index));
    }
  } catch (std::exception &e) {
    parser.error.insert(e.what());
    return;
  }
  if (parser.valid()) {
    parser.value = std::move(p);
  }
}

}  // namespace

// OccPosition

jsonParser &to_json(
    occ_events::OccPosition const &pos, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system) {
  json.put_obj();
  if (pos.is_in_reservoir) {
    json["is_in_reservoir"] = true;
    json["occupant_index"] = pos.occupant_index;
    if (system.has_value()) {
      json["chemical_name"] = system->get().get_chemical_name(pos);
    }
  } else {
    json["coordinate"] = pos.integral_site_coordinate;
    json["occupant_index"] = pos.occupant_index;
    if (!pos.is_atom) {
      json["molecule"] = true;
    } else {
      json["atom_position_index"] = pos.atom_position_index;
    }

    if (system.has_value()) {
      std::string chemical_name = system->get().get_chemical_name(pos);
      json["chemical_name"] = system->get().get_chemical_name(pos);
      std::string orientation_name = system->get().get_orientation_name(pos);
      if (chemical_name != orientation_name) {
        json["orientation_name"] = orientation_name;
      }
      if (!pos.is_atom) {
        Index b = pos.integral_site_coordinate.sublattice();
        Index mol_size = system->get()
                             .atom_position_to_name_index[b][pos.occupant_index]
                             .size();
        if (mol_size > 1) {
          json["atom_name"] = system->get().get_atom_name(pos);
        }
      }
    }
  }
  return json;
}

void from_json(occ_events::OccPosition &pos, jsonParser const &json,
               occ_events::OccSystem const &system) {
  pos = jsonConstructor<occ_events::OccPosition>::from_json(json, system);
}

occ_events::OccPosition jsonConstructor<occ_events::OccPosition>::from_json(
    jsonParser const &json, occ_events::OccSystem const &system) {
  InputParser<occ_events::OccPosition> parser{json, system};
  std::stringstream ss;
  ss << "Error: Invalid OccPosition JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

void parse(InputParser<occ_events::OccPosition> &parser,
           occ_events::OccSystem const &system) {
  // parse "is_in_reservoir"
  bool is_in_reservoir = false;
  parser.optional(is_in_reservoir, "is_in_reservoir");

  if (is_in_reservoir) {
    parse_reservoir_position(parser, system);
  } else {
    parse_non_reservoir_position(parser, system);
  }
}

// OccTrajectory

jsonParser &to_json(
    occ_events::OccTrajectory const &traj, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system) {
  json.put_array();
  for (auto const &pos : traj.position) {
    json.push_back(pos, system);
  }
  return json;
}

void from_json(occ_events::OccTrajectory &traj, jsonParser const &json,
               occ_events::OccSystem const &system) {
  traj = jsonConstructor<occ_events::OccTrajectory>::from_json(json, system);
}

occ_events::OccTrajectory jsonConstructor<occ_events::OccTrajectory>::from_json(
    jsonParser const &json, occ_events::OccSystem const &system) {
  InputParser<occ_events::OccTrajectory> parser{json, system};
  std::stringstream ss;
  ss << "Error: Invalid OccTrajectory JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

void parse(InputParser<occ_events::OccTrajectory> &parser,
           occ_events::OccSystem const &system) {
  if (!parser.self.is_array()) {
    parser.error.insert("Error: expected JSON array for OccTrajectory");
    return;
  }

  parser.value = notstd::make_unique<occ_events::OccTrajectory>();
  occ_events::OccTrajectory &traj = *parser.value;
  for (Index i = 0; i < parser.self.size(); ++i) {
    fs::path path{std::to_string(i)};
    auto shared_subparser =
        parser.subparse<occ_events::OccPosition>(path, system);
    if (!shared_subparser->valid()) {
      parser.value.reset();
      return;
    }
    traj.position.push_back(*shared_subparser->value);
  }
}

// OccEvent

jsonParser &to_json(
    occ_events::OccEvent const &event, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system,
    occ_events::OccEventOutputOptions const &options) {
  json["trajectories"].put_array();
  for (auto const &traj : event) {
    json["trajectories"].push_back(traj, system);
  }

  clust::IntegralCluster clust = make_cluster(event);

  if (options.include_cluster_occupation) {
    try {
      auto cluster_occupation = make_cluster_occupation(event);
      json["occupation"]["index"]["initial"] = cluster_occupation.second[0];
      json["occupation"]["index"]["final"] = cluster_occupation.second[1];

      if (system.has_value()) {
        to_json(clust, json["cluster"], *system->get().prim);

        auto _names = [&](std::vector<int> const &occ) {
          std::vector<std::string> n;
          Index i = 0;
          for (auto const &site : clust) {
            n.push_back(system->get().get_orientation_name(site, occ[i]));
            ++i;
          }
          return n;
        };
        json["occupation"]["name"]["initial"] =
            _names(cluster_occupation.second[0]);
        json["occupation"]["name"]["final"] =
            _names(cluster_occupation.second[1]);
      }
    } catch (std::exception &e) {
      if (system.has_value()) {
        clust::IntegralCluster clust = make_cluster(event);
        to_json(clust, json["cluster"], *system->get().prim);
      }
      json["occupation"]["error"] = "inconsistent occupanacies";
    }
  } else if (options.include_cluster) {
    if (system.has_value()) {
      to_json(clust, json["cluster"], *system->get().prim);
    }
  }
  if (options.include_event_invariants) {
    if (system.has_value()) {
      occ_events::OccEventInvariants invariants(event, system->get());
      to_json(invariants, json["event_invariants"]);
    }
  }
  return json;
}

/// \brief OccEvent orbit printing
jsonParser &to_json(
    std::set<occ_events::OccEvent> const &orbit, jsonParser &json,
    std::optional<std::reference_wrapper<occ_events::OccSystem const>> system,
    std::shared_ptr<occ_events::SymGroup const> const &factor_group,
    std::vector<occ_events::OccEventRep> const &occevent_symgroup_rep,
    occ_events::OccEventOutputOptions const &options) {
  // --- "mult" ---
  json["mult"] = orbit.size();
  if (!orbit.size()) {
    return json;
  }

  // --- "prototype" ---
  to_json(*orbit.begin(), json["prototype"], system, options);

  // --- "elements" ---
  if (options.include_elements) {
    json["elements"].put_array();
    for (auto const &event : orbit) {
      jsonParser tjson;
      to_json(event, tjson, system, options);
      json["elements"].push_back(tjson);
    }
  }

  // --- lambda to print symop descriptions ---
  auto _put_desc = [&](std::vector<Index> const &factor_group_indices,
                       jsonParser &j) {
    j.put_array();
    for (Index i : factor_group_indices) {
      xtal::SymInfo syminfo{factor_group->element[i],
                            system->get().prim->lattice()};
      j.push_back(to_brief_unicode(syminfo, options.sym_info_options));
    }
    j.set_multiline_array();
  };

  // --- invariant_group ---
  if (options.include_invariant_group && system.has_value()) {
    std::vector<std::shared_ptr<occ_events::SymGroup const>> g =
        make_occevent_groups(orbit, factor_group,
                             system->get().prim->lattice().lat_column_mat(),
                             occevent_symgroup_rep);

    to_json(g[0]->head_group_index, json["prototype"]["invariant_group"]);
    _put_desc(g[0]->head_group_index,
              json["prototype"]["invariant_group_descriptions"]);

    if (options.include_elements) {
      for (Index i = 0; i < orbit.size(); ++i) {
        to_json(g[i]->head_group_index, json["elements"][i]["invariant_group"]);
        _put_desc(g[i]->head_group_index,
                  json["elements"][i]["invariant_group_descriptions"]);
      }
    }
  }

  // --- equivalence_map ---
  if (options.include_equivalence_map && system.has_value()) {
    std::vector<std::vector<Index>> eq_map = group::make_equivalence_map(
        orbit, occevent_symgroup_rep.begin(), occevent_symgroup_rep.end(),
        occ_events::prim_periodic_occevent_copy_apply);

    to_json(eq_map[0], json["prototype"]["equivalence_map"]);
    _put_desc(eq_map[0], json["prototype"]["equivalence_map_descriptions"]);

    if (options.include_elements) {
      for (Index i = 0; i < orbit.size(); ++i) {
        to_json(eq_map[i], json["elements"][i]["equivalence_map"]);
        _put_desc(eq_map[i],
                  json["elements"][i]["equivalence_map_descriptions"]);
      }
    }
  }
  return json;
}

void from_json(occ_events::OccEvent &event, jsonParser const &json,
               occ_events::OccSystem const &system) {
  event = jsonConstructor<occ_events::OccEvent>::from_json(json, system);
}

occ_events::OccEvent jsonConstructor<occ_events::OccEvent>::from_json(
    jsonParser const &json, occ_events::OccSystem const &system) {
  InputParser<occ_events::OccEvent> parser{json, system};
  std::stringstream ss;
  ss << "Error: Invalid OccEvent JSON object";
  report_and_throw_if_invalid(parser, err_log(), std::runtime_error{ss.str()});
  return *parser.value;
}

void parse(InputParser<occ_events::OccEvent> &parser,
           occ_events::OccSystem const &system) {
  if (!parser.self.is_obj()) {
    parser.error.insert(
        "Error: expected JSON object for OccEvent 'trajectories'");
    return;
  }

  if (!parser.self.contains("trajectories")) {
    parser.error.insert("Error: expected JSON array for OccEvent");
    return;
  }

  if (!parser.self["trajectories"].is_array()) {
    parser.error.insert(
        "Error: expected JSON array for OccEvent 'trajectories'");
    return;
  }

  parser.value = notstd::make_unique<occ_events::OccEvent>();
  occ_events::OccEvent &event = *parser.value;
  for (Index i = 0; i < parser.self["trajectories"].size(); ++i) {
    auto shared_subparser = parser.subparse<occ_events::OccTrajectory>(
        fs::path("trajectories") / std::to_string(i), system);
    if (!shared_subparser->valid()) {
      parser.value.reset();
      return;
    }
    event.elements().push_back(*shared_subparser->value);
  }
}

// OccEventInvariants

jsonParser &to_json(occ_events::OccEventInvariants const &invariants,
                    jsonParser &json) {
  json["distances"] = invariants.distances();
  json["molecule_count"].put_array();
  for (auto const &count : invariants.molecule_count()) {
    jsonParser tjson;
    to_json(count, tjson, jsonParser::as_array());
    json["molecule_count"].push_back(tjson);
  }
  json["size"] = invariants.size();
  return json;
}

}  // namespace CASM
