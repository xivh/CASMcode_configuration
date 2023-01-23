#include "casm/configuration/enumeration/MakeOccEventStructures.hh"

#include "casm/configuration/Configuration.hh"
#include "casm/configuration/clusterography/IntegralCluster.hh"
#include "casm/configuration/occ_events/OccEvent.hh"
#include "casm/configuration/occ_events/OccSystem.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/StrainConverter.hh"

// debug:
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"

namespace CASM {
namespace config {

/// \brief MakeOccEventStructures constructor
///
/// \param configuration Configuration specifying the occupation
///     on all non-event sites
/// \param occ_event The event
/// \param system System used to specify order and indices
/// \param skip_event_occupants If true, do not include occupants
///     involved in the event in the output. Default=false.
///
MakeOccEventStructures::MakeOccEventStructures(
    Configuration const &configuration, occ_events::OccEvent const &occ_event,
    std::shared_ptr<occ_events::OccSystem const> const &system,
    bool skip_event_occupants) {
  // references
  auto const &supercell = *configuration.supercell;
  auto const &prim = *supercell.prim;
  auto const &converter = supercell.unitcellcoord_index_converter;
  auto const &dof_values = configuration.dof_values;
  auto const &occupation = dof_values.occupation;
  auto const &local_dof_values = dof_values.local_dof_values;
  auto const &global_dof_values = dof_values.global_dof_values;

  // validate no multi-atom molecules
  std::vector<xtal::Molecule> molecule_list =
      xtal::struc_molecule(*prim.basicstructure);
  for (auto const &mol : molecule_list) {
    if (mol.size() > 1) {
      throw std::runtime_error(
          "Error in MakeOccEventStructures:  multi-atom molecules are not "
          "allowed");
    }
  }

  // validate no intermediate event states
  for (auto const &traj : occ_event) {
    if (traj.position.size() != 2) {
      throw std::runtime_error(
          "Error in MakeOccEventStructures: a trajectory does not have 2 "
          "positions");
    }
  }

  // set m_ideal_lat_column_mat
  m_ideal_lat_column_mat =
      supercell.superlattice.superlattice().lat_column_mat();

  // set m_F
  if (has_strain_dof(*prim.basicstructure)) {
    if (global_dof_values.size() > 1) {
      throw std::runtime_error(
          "Error in MakeOccEventStructures: only occupation and strain DoF are "
          "supported");
    }
    DoFKey dof_key = get_strain_dof_key(*prim.basicstructure);
    xtal::StrainConverter DoFstrain_converter(
        dof_key, prim.global_dof_info.at(dof_key).basis());
    m_F = DoFstrain_converter.to_F(global_dof_values.at(dof_key));
  } else {
    if (global_dof_values.size() > 0) {
      throw std::runtime_error(
          "Error in MakeOccEventStructures: only occupation strain DoF are "
          "supported");
    }
    m_F = Eigen::Matrix3d::Identity();
  }

  // validate no local DoF
  if (local_dof_values.size() > 0) {
    throw std::runtime_error(
        "Error in MakeOccEventStructures: only occupation strain DoF are "
        "supported");
  }

  // here we set m_atom_names, m_coords, and m_disp:
  auto cluster_occupation = make_cluster_occupation(occ_event);
  std::set<Index> event_sites =
      to_index_set(cluster_occupation.first, converter);
  for (Index chemical_index = 0;
       chemical_index < system->chemical_name_list.size(); ++chemical_index) {
    // skip vacancies
    if (system->is_vacancy_list[chemical_index]) {
      continue;
    }
    std::string chemical_name = system->chemical_name_list[chemical_index];

    if (skip_event_occupants == false) {
      // add event occupants if matching current chemical_name
      for (occ_events::OccTrajectory const &traj : occ_event) {
        occ_events::OccPosition const &pos_init = traj.position[0];
        if (system->get_chemical_index(pos_init) == chemical_index) {
          occ_events::OccPosition const &pos_final = traj.position[1];
          Eigen::Vector3d coord_init =
              system->get_cartesian_coordinate(pos_init);
          Eigen::Vector3d coord_final =
              system->get_cartesian_coordinate(pos_final);
          m_atom_names.push_back(chemical_name);
          m_coords_init.push_back(coord_init);
          m_event_disp.push_back(coord_final - coord_init);
        }
      }
    }

    // add configuration occupations,
    // if matching current chemical_name and not on an event site
    for (Index l = 0; l < converter.total_sites(); ++l) {
      if (event_sites.count(l)) {
        continue;
      }
      occ_events::OccPosition pos =
          system->make_molecule_position(converter(l), occupation(l));
      if (system->get_chemical_index(pos) == chemical_index) {
        m_atom_names.push_back(chemical_name);
        m_coords_init.push_back(system->get_cartesian_coordinate(pos));
        m_event_disp.push_back(Eigen::Vector3d::Zero());
      }
    }
  }
}

/// \brief Construct interpolated structure
///
/// \param interpolation_factor Interpolation factor, where 0.0 results
///     in the initial configuration and 1.0 results in the final
///     configuration. Other values linearly interpolate the
///     displacements of the atoms involved in the event.
///
/// \return structure The interpolated structure
///
xtal::SimpleStructure MakeOccEventStructures::operator()(
    double interpolation_factor) const {
  xtal::SimpleStructure structure;
  structure.lat_column_mat = m_ideal_lat_column_mat;
  structure.atom_info.resize(m_atom_names.size());
  structure.atom_info.names = m_atom_names;
  double f = interpolation_factor;
  for (Index i = 0; i < m_atom_names.size(); ++i) {
    structure.atom_info.coords.col(i) = m_coords_init[i] + f * m_event_disp[i];
  }
  structure.deform_coords(m_F);
  return structure;
}

}  // namespace config
}  // namespace CASM
