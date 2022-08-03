#include "casm/configuration/sym_info/unitcellcoord_sym_info.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/UnitCellCoordRep.hh"

namespace CASM {
namespace sym_info {

UnitCellCoordSymGroupRep make_unitcellcoord_symgroup_rep(
    std::vector<xtal::SymOp> const &group_elements,
    xtal::BasicStructure const &basicstructure) {
  // sitemap[i]->std::vector<UnitCellCoord>:
  // - the vector of sites that i-th op maps prim basis sites onto
  // - used for constructing the site-related symreps
  std::vector<std::vector<UnitCellCoord>> sitemap;
  for (SymOp const &op : group_elements) {
    sitemap.push_back(xtal::symop_site_map(op, basicstructure));
  }

  // Construct (BasisPermutationSymGroupRep) unitcellcoord_symgroup_rep
  Index op_index = 0;
  UnitCellCoordSymGroupRep unitcellcoord_symgroup_rep;
  for (SymOp const &op : group_elements) {
    UnitCellCoordRep rep;
    rep.point_matrix = lround(cart2frac(op.matrix, basicstructure.lattice()));
    for (UnitCellCoord const &site : sitemap[op_index]) {
      rep.sublattice_index.push_back(site.sublattice());
      rep.unitcell_indices.push_back(site.unitcell());
    }
    unitcellcoord_symgroup_rep.push_back(rep);
    ++op_index;
  }
  return unitcellcoord_symgroup_rep;
}

}  // namespace sym_info
}  // namespace CASM
