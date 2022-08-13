#include "casm/configuration/occ_events/io/json/OccSystem_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/configuration/occ_events/OccSystem.hh"

namespace CASM {

jsonParser &to_json(occ_events::OccSystem const &system, jsonParser &json) {
  to_json(system.chemical_name_list, json["chemical_name_list"]);
  to_json(system.is_indivisible_chemical_list,
          json["is_indivisible_chemical_list"]);
  to_json(system.vacancy_name_list, json["vacancy_name_list"]);
  to_json(system.is_vacancy_list, json["is_vacancy_list"]);
  to_json(system.atom_name_list, json["atom_name_list"]);
  to_json(system.orientation_name_list, json["orientation_name_list"]);
  to_json(system.atom_position_to_name_index,
          json["atom_position_to_name_index"]);
  to_json(system.occupant_to_chemical_index,
          json["occupant_to_chemical_index"]);
  to_json(system.occupant_to_orientation_index,
          json["occupant_to_orientation_index"]);
  to_json(system.orientation_to_occupant_index,
          json["orientation_to_occupant_index"]);
  return json;
}

}  // namespace CASM
