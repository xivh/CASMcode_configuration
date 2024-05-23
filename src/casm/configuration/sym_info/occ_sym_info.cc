#include "casm/configuration/sym_info/occ_sym_info.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/OccupantDoFIsEquivalent.hh"
#include "casm/crystallography/UnitCellCoord.hh"

namespace CASM {
namespace sym_info {

OccSymInfo::OccSymInfo(std::vector<xtal::SymOp> const &group_elements,
                       xtal::BasicStructure const &basicstructure) {
  // sitemap[i]->std::vector<UnitCellCoord>:
  // - the vector of sites that i-th op maps prim basis sites onto
  // - used for constructing the site-related symreps
  std::vector<std::vector<UnitCellCoord>> sitemap;
  for (SymOp const &op : group_elements) {
    sitemap.push_back(xtal::symop_site_map(op, basicstructure));
  }

  auto const &basis = basicstructure.basis();

  // Construct (OccSymGroupRep) occ_symgroup_rep, atom_position_rep,
  //   & set (bool) has_occupation_dofs & set (bool) has_aniso_occs
  this->has_occupation_dofs = false;
  for (Index b = 0; b < basis.size(); ++b) {
    if (basis[b].occupant_dof().size() > 1) {
      this->has_occupation_dofs = true;
    }
  }

  this->has_aniso_occs = false;
  Index op_index = 0;
  for (SymOp const &op : group_elements) {
    // the occ_rep is a vector of permutation (one per prim site), that are used
    // to generate the new occupant indices
    OccSymOpRep occ_rep;

    // atom_position_rep[sublattice_index_before][occupant_index_before] ->
    // Permutation
    AtomPositionSymOpRep atom_position_rep;

    for (Index b = 0; b < basis.size(); ++b) {
      // copy_aply(symop,dofref_from) = P.permute(dofref_to);
      auto const &dofref_to =
          basis[sitemap[op_index][b].sublattice()].occupant_dof();
      auto const &dofref_from = basis[b].occupant_dof();

      xtal::OccupantDoFIsEquivalent eq(dofref_from);

      if (eq(op, dofref_to)) {
        if (!eq.perm().is_identity()) {
          this->has_aniso_occs = true;
        }
        occ_rep.push_back(eq.perm().perm_array());

        std::vector<Permutation> tmp;
        for (Index i = 0; i < dofref_from.size(); ++i) {
          tmp.push_back(eq.atom_position_perm()[i].perm_array());
        }
        atom_position_rep.push_back(tmp);
      } else {
        throw std::runtime_error(
            "In CASM::config::Prim constructor: Sites originally "
            "identified as equivalent cannot be mapped by symmetry.");
      }
    }
    this->occ_symgroup_rep.push_back(occ_rep);
    this->atom_position_symgroup_rep.push_back(atom_position_rep);

    ++op_index;
  }
}

}  // namespace sym_info
}  // namespace CASM
