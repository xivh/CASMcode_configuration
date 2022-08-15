#include "casm/configuration/copy_configuration.hh"

#include "casm/configuration/ConfigIsEquivalent.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/SupercellSymOp.hh"
#include "casm/configuration/canonical_form.hh"

namespace CASM {
namespace config {

namespace {

SupercellSymOp find_translation(Configuration const &configuration) {
  auto begin = SupercellSymOp::translation_begin(configuration.supercell);
  auto end = SupercellSymOp::translation_end(configuration.supercell);
  if (++begin == end) {
    return end;
  }
  ConfigIsEquivalent f(configuration);
  return std::find_if(begin, end, f);
}

}  // namespace

/// \brief Copy configuration DoF values into a supercell
///
/// \param motif The initial configuration
/// \param supercell The Supercell of the new configuration
/// \param origin The UnitCell indicating which unit cell in the
///        initial configuration is the origin in new configuration
///
/// Notes:
/// - This method assumes the motif forms an infinite crystal and copies site
///   DoF values that lie inside `supercell` directory into a new configuration.
///
Configuration copy_configuration(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell, UnitCell const &origin) {
  if (supercell->prim != motif.supercell->prim) {
    throw std::runtime_error(
        "Error in CASM::config::copy_configuration: prim mismatch.");
  }

  Index supercell_total_sites =
      supercell->unitcellcoord_index_converter.total_sites();

  // construct new_config
  Configuration new_config{supercell};

  // copy global DoF values
  new_config.dof_values.global_dof_values = motif.dof_values.global_dof_values;

  // copy occupation values
  for (Index i = 0; i < supercell_total_sites; i++) {
    // unitcellcoord of site i in new_config
    UnitCellCoord unitcellcoord =
        new_config.supercell->unitcellcoord_index_converter(i);

    // equivalent site in configuration
    Index motif_site_index =
        motif.supercell->unitcellcoord_index_converter(unitcellcoord + origin);

    // copy occupation value
    new_config.dof_values.occupation(i) =
        motif.dof_values.occupation(motif_site_index);
  }

  // copy local DoF values
  for (auto const &pair : motif.dof_values.local_dof_values) {
    std::string name = pair.first;
    Eigen::MatrixXd const &M_config = pair.second;
    Eigen::MatrixXd &M_sub = new_config.dof_values.local_dof_values.at(name);
    for (Index i = 0; i < supercell_total_sites; i++) {
      // unitcellcoord of site i in new_config
      UnitCellCoord unitcellcoord =
          new_config.supercell->unitcellcoord_index_converter(i);

      // equivalent site in superconfig
      Index motif_site_index = motif.supercell->unitcellcoord_index_converter(
          unitcellcoord + origin);

      // copy dof from superconfig to this:
      M_sub.col(i) = M_config.col(motif_site_index);
    }
  }

  return new_config;
}

/// \brief Copy transformed configuration DoF values into a supercell
///
/// \param prim_factor_group_index Index of prim factor group operation which
///     transforms the initial configuration
/// \param translation Lattice translation applied after the prim factor group
///     operation
/// \param motif The initial configuration
/// \param supercell The Supercell of the new configuration
/// \param origin The UnitCell indicating which unit cell in the
///        transformed configuration is the origin in new configuration
///
/// Copies DoF values as if `motif` is transformed by the prim factor group
/// operation with index `factor_group_index`, then translated by
/// `translation`, then copied starting from `origin`. In other words, sites
/// map according to:
///     new_config_unitcellcoord + origin = fg * motif_unitcellcoord + trans
///
Configuration copy_configuration(
    Index prim_factor_group_index, UnitCell translation,
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell, UnitCell const &origin) {
  if (supercell->prim != motif.supercell->prim) {
    throw std::runtime_error(
        "Error in CASM::config::copy_configuration (and transform): prim "
        "mismatch.");
  }

  auto const &prim = supercell->prim;
  PrimSymInfo const &prim_sym_info = prim->sym_info;
  auto const &unitcellcoord_rep = prim_sym_info.unitcellcoord_symgroup_rep;
  Index inverse_prim_factor_group_index =
      prim->sym_info.factor_group->inverse_index[prim_factor_group_index];

  Index supercell_total_sites =
      supercell->unitcellcoord_index_converter.total_sites();

  // construct new_config
  Configuration new_config{supercell};

  // copy transformed global DoF values
  for (auto const &pair : motif.dof_values.global_dof_values) {
    std::string name = pair.first;
    Eigen::VectorXd const &V_motif = pair.second;
    Eigen::VectorXd &V_new = new_config.dof_values.global_dof_values.at(name);
    auto const &global_rep = prim_sym_info.global_dof_symgroup_rep.at(name);

    Eigen::MatrixXd const &matrix_rep = global_rep[prim_factor_group_index];
    V_new = matrix_rep * V_motif;
  }

  // unitcellcoord + origin = fg * motif_unitcellcoord + trans
  // motif_unitcellcoord = fg_inverse * (unitcellcoord  + origin - trans)

  // copy transformed occupation values
  for (Index i = 0; i < supercell_total_sites; i++) {
    // unitcellcoord of site i in new_config
    UnitCellCoord unitcellcoord =
        new_config.supercell->unitcellcoord_index_converter(i);

    // motif_unitcellcoord, the site which transforms to site i in new_config
    UnitCellCoord motif_unitcellcoord =
        copy_apply(unitcellcoord_rep[inverse_prim_factor_group_index],
                   (unitcellcoord + origin - translation));

    // equivalent site in configuration
    Index motif_site_index =
        motif.supercell->unitcellcoord_index_converter(motif_unitcellcoord);

    // occupation value transformation (accounts for aniostropic occupants)
    Index b = motif_unitcellcoord.sublattice();
    sym_info::Permutation const &perm_rep =
        prim_sym_info.occ_symgroup_rep[prim_factor_group_index][b];

    // copy transformed occupation value
    new_config.dof_values.occupation(i) =
        perm_rep[motif.dof_values.occupation(motif_site_index)];
  }

  // copy local DoF values
  for (auto const &pair : motif.dof_values.local_dof_values) {
    std::string name = pair.first;
    Eigen::MatrixXd const &M_motif = pair.second;
    Eigen::MatrixXd &M_new = new_config.dof_values.local_dof_values.at(name);
    auto const &local_rep = prim_sym_info.local_dof_symgroup_rep.at(name);
    for (Index i = 0; i < supercell_total_sites; i++) {
      // unitcellcoord of site i in new_config
      UnitCellCoord unitcellcoord =
          new_config.supercell->unitcellcoord_index_converter(i);

      // motif_unitcellcoord, the site which transforms to site i in new_config
      UnitCellCoord motif_unitcellcoord =
          copy_apply(unitcellcoord_rep[inverse_prim_factor_group_index],
                     (unitcellcoord + origin - translation));

      // equivalent site in configuration
      Index motif_site_index =
          motif.supercell->unitcellcoord_index_converter(motif_unitcellcoord);

      // local DoF value transformation
      Index b = motif_unitcellcoord.sublattice();
      Eigen::MatrixXd const &matrix_rep = local_rep[prim_factor_group_index][b];

      // copy dof from superconfig to this:
      M_new.col(i) = matrix_rep * M_motif.col(motif_site_index);
    }
  }

  return new_config;
}

/// \brief Make all equivalent configurations with respect to the prim factor
/// group
///     that fill a supercell
///
/// \param motif The motif configuration
/// \param supercell The supercell to fill
///
/// \returns All configurations equivalent with respect to the prim factor
///     group which fit in the given supercell.
std::vector<Configuration> make_all_super_configurations(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell) {
  Configuration prim_motif = make_primitive(motif);
  xtal::Lattice prim_motif_lattice =
      prim_motif.supercell->superlattice.superlattice();
  xtal::Lattice supercell_lattice = supercell->superlattice.superlattice();
  Prim const &prim = *motif.supercell->prim;
  SymGroup const &prim_fg = *prim.sym_info.factor_group;
  SymGroup const &supercell_fg = *supercell->sym_info.factor_group;
  SymGroup const &prim_motif_supercell_fg =
      *prim_motif.supercell->sym_info.factor_group;
  double xtal_tol = prim.basicstructure->lattice().tol();

  // - Want to find the unique ways to fill supercell with prim_motif.
  // - Will be doing prim_fg_op * prim_motif, but only if
  //   supercell_lattice is a supercell of prim_fg_op*prim_motif_lattice.
  // - To only keep the unique ways of re-orienting, given a prim_fg_op,
  //   keep only one of the set generated by:
  //
  //        supercell_fg_op * prim_fg_op * prim_motif_supercell_fg_op.
  //
  // - So keep prim_fg_op if it is the minimum of all the combined ops

  auto generates_unique_orientation = [&](Index prim_fg_op) {
    for (Index supercell_fg_op : supercell_fg.head_group_index) {
      for (Index prim_motif_scel_fg_op :
           prim_motif_supercell_fg.head_group_index) {
        Index combined_op = prim_fg.mult(
            supercell_fg_op, prim_fg.mult(prim_fg_op, prim_motif_scel_fg_op));
        if (combined_op < prim_fg_op) {
          return false;
        }
      }
    }
    return true;
  };

  std::set<Index> unique_generating_prim_fg_op;
  for (Index i = 0; i < prim_fg.element.size(); ++i) {
    if (generates_unique_orientation(i)) {
      unique_generating_prim_fg_op.insert(i);
    }
  }

  std::vector<Configuration> all;
  UnitCell trans(0, 0, 0);
  UnitCell origin(0, 0, 0);
  SupercellSymOp begin = SupercellSymOp::begin(supercell);
  SupercellSymOp end = SupercellSymOp::end(supercell);

  // Loop over unique generating ops
  for (Index prim_fg_op : unique_generating_prim_fg_op) {
    // If prim_fg_op * prim_motif doesn't fill supercell, skip
    auto test_lattice =
        sym::copy_apply(prim_fg.element[prim_fg_op], prim_motif_lattice);
    if (!is_superlattice(supercell_lattice, test_lattice, xtal_tol).first) {
      continue;
    }

    // Apply op to fill supercell and make all equivalents
    Configuration tmp =
        copy_configuration(prim_fg_op, trans, prim_motif, supercell, origin);
    std::vector<Configuration> equiv = make_equivalents(tmp, begin, end);
    all.insert(std::end(all), std::begin(equiv), std::end(equiv));
  }
  return all;
}

/// \brief Make all equivalent configurations with respect to the prim factor
/// group
///     that fill a supercell
///
/// \param motif The motif configuration
/// \param supercell The supercell to fill
///
/// \returns All configurations equivalent with respect to the prim factor
///     group which fit in the given supercell.
std::set<Configuration> make_all_super_configurations_check(
    Configuration const &motif,
    std::shared_ptr<Supercell const> const &supercell) {
  Configuration prim_motif = make_primitive(motif);
  xtal::Lattice prim_motif_lattice =
      prim_motif.supercell->superlattice.superlattice();
  xtal::Lattice supercell_lattice = supercell->superlattice.superlattice();
  Prim const &prim = *motif.supercell->prim;
  SymGroup const &prim_fg = *prim.sym_info.factor_group;
  double xtal_tol = prim.basicstructure->lattice().tol();

  std::set<Configuration> all;
  UnitCell trans(0, 0, 0);
  UnitCell origin(0, 0, 0);
  SupercellSymOp begin = SupercellSymOp::begin(supercell);
  SupercellSymOp end = SupercellSymOp::end(supercell);

  // Loop over prim factor group ops
  for (Index prim_fg_op = 0; prim_fg_op < prim_fg.element.size();
       ++prim_fg_op) {
    // If prim_fg_op * prim_motif doesn't fill supercell, skip
    auto test_lattice =
        sym::copy_apply(prim_fg.element[prim_fg_op], prim_motif_lattice);
    if (!is_superlattice(supercell_lattice, test_lattice, xtal_tol).first) {
      continue;
    }

    // Apply op to fill supercell, then make all equivalents
    Configuration tmp =
        copy_configuration(prim_fg_op, trans, prim_motif, supercell, origin);
    if (!all.count(tmp)) {
      for (auto it = begin; it != end; ++it) {
        all.emplace(copy_apply(*it, tmp));
      }
    }
  }
  return all;
}

/// \brief Return true if no translations within the supercell result in the
///     same configuration
bool is_primitive(Configuration const &configuration) {
  auto end = SupercellSymOp::translation_end(configuration.supercell);
  return find_translation(configuration) == end;
}

/// \brief Return the primitive configuration
///
/// Notes:
/// - Does not apply any symmetry operations
/// - Use `make_in_canonical_supercell` aftwards to obtain the primitive
///   canonical configuration in the canonical supercell.
Configuration make_primitive(Configuration const &configuration) {
  Configuration tconfig{configuration};
  auto prim = tconfig.supercell->prim;
  double xtal_tol = prim->basicstructure->lattice().tol();

  // check if config is primitive, and if not, obtain a translation that maps
  // the config on itself
  while (true) {
    SupercellSymOp result = find_translation(tconfig);

    if (result == SupercellSymOp::translation_end(tconfig.supercell)) {
      break;
    }

    // replace one of the lattice vectors with the translation
    Lattice new_lat =
        xtal::replace_vector(tconfig.supercell->superlattice.superlattice(),
                             result.to_symop().translation, xtal_tol)
            .make_right_handed()
            .reduced_cell();

    // create a sub configuration in the new supercell
    tconfig =
        copy_configuration(tconfig, std::make_shared<Supercell>(prim, new_lat));
  }

  return tconfig;
}

/// \brief Return the canonical configuration in the canonical supercell
///
/// Notes:
/// - Applies symmetry operations if necessary to copy the configuration into
///   the canonical supercell
Configuration make_in_canonical_supercell(Configuration const &configuration) {
  if (is_canonical(*configuration.supercell)) {
    return make_canonical_form(configuration,
                               SupercellSymOp::begin(configuration.supercell),
                               SupercellSymOp::end(configuration.supercell));
  }

  std::shared_ptr<Supercell const> canonical_supercell =
      make_canonical_form(*configuration.supercell);
  Lattice const &canon_scel_lattice =
      canonical_supercell->superlattice.superlattice();
  Lattice const &scel_lattice =
      configuration.supercell->superlattice.superlattice();
  auto const &prim_fg = *configuration.supercell->prim->sym_info.factor_group;
  auto begin = prim_fg.element.begin();
  auto end = prim_fg.element.end();
  auto res = xtal::is_equivalent_superlattice(canon_scel_lattice, scel_lattice,
                                              begin, end, scel_lattice.tol());
  Index prim_factor_group_index = std::distance(begin, res.first);

  Configuration config_in_canonical_supercell = copy_configuration(
      prim_factor_group_index, {0, 0, 0}, configuration, canonical_supercell);

  return make_canonical_form(
      config_in_canonical_supercell,
      SupercellSymOp::begin(config_in_canonical_supercell.supercell),
      SupercellSymOp::end(config_in_canonical_supercell.supercell));
}

}  // namespace config
}  // namespace CASM
