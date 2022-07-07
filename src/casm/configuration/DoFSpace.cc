#include "casm/configuration/DoFSpace.hh"

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/clexulator/ConfigDoFValuesTools_impl.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
namespace config {

namespace DoFSpace_impl {

// defined in DoFSpaceAxisInfo.cc
void throw_if_missing_local_dof_requirements(
    DoFKey const &dof_key,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites);

}  // namespace DoFSpace_impl

/// \brief DoFSpace
///
/// The DoFSpace class is used to specify a particular degree of freedom space.
/// DoFSpace have:
///   - dof_key (DoFKey): A string indicating which DoF type (e.g., "disp",
///   "Hstrain", "occ")
///   - prim (std::shared_ptr<xtal::BasicStructure const>): The prim structure
///   - transformation_matrix_to_super (std::optional<Eigen::Matrix3l>):
///   Specifies the supercell for a local DoF space. Has value for local DoF
///   only.
///   - sites (std::optional<std::set<Index>>): The sites included in a local
///   DoF space. Has value for local DoF only.
///   - basis (Eigen::MatrixXd): Allows specifying a subspace of the space
///   determined from dof_key, and for local DoF, transformation_matrix_to_super
///   and sites. The rows of `basis` correspond to prim DoF basis axes, so the
///   following relations apply:
///
///       standard_dof_values = dof_set.basis() * prim_dof_values
///       prim_dof_values = dof_space.basis() * dof_space_coordinate
///
///    Examples:
///
///       For dof_key=="disp", and sites.value().size()==4 (any
///       transformation_matrix_to_super):
///         The default DoF space has dimension 12 corresponding to (x,y,z) for
///         each of the 4 sites. Then, basis is a matrix with 12 rows and 12 or
///         fewer columns.
///
///       For dof_key=="occ", and sites().value().size()==4 (any
///       transformation_matrix_to_super), with the particular sites selected
///       having allowed occupations ["A", "B"], ["A", "B", "C"], ["A, B"], and
///       ["A", "B", "C"] :
///         The default DoF space has dimension 10 = 2 + 3 + 2 + 3.
///         Then, basis is a matrix with 10 rows and 10 or fewer columns.
///
///       For dof_key=="GLstrain":
///         The default DoF space has dimension 6 for the six independent
///         GLstrain components. Then, basis is a matrix with 6 rows and 6 or
///         fewer columns.
///
///     By default, basis is the full DoF space (identity matrix with dimensions
///     matching the full space for the particular combination of config_region
///     and dof_key).
///

/// \brief Constructor
///
/// \param dof_key (DoFKey): A string indicating which DoF type (e.g., "disp",
///   "Hstrain", "occ")
/// \param prim (std::shared_ptr<xtal::BasicStructure const>): The prim
/// structure \param transformation_matrix_to_super
/// (std::optional<Eigen::Matrix3l>):
///   Specifies the supercell for a local DoF space. Ignored for global DoF.
/// \param sites (std::optional<std::set<Index>>): The sites included in a local
///   DoF space. Ignored for global DoF. If dof_key specifies a local DoF and
///   this does not have a value, all sites in the supercell are incldued.
/// \param basis (Eigen::MatrixXd): Allows specifying a subspace of the space
///   determined from dof_key, and for local DoF, transformation_matrix_to_super
///   and sites. The rows of `basis` correspond to prim DoF basis axes, see
///   class documentation for relations apply. If this does not have a value,
///   the standard basis (identify matrix of appropriate dimension) is used.
///
/// \seealso make_dof_space
DoFSpace::DoFSpace(
    DoFKey const &_dof_key,
    std::shared_ptr<xtal::BasicStructure const> const &_prim,
    std::optional<Eigen::Matrix3l> _transformation_matrix_to_super,
    std::optional<std::set<Index>> _sites,
    std::optional<Eigen::MatrixXd> _basis)
    : DoFSpace(make_dof_space(_dof_key, _prim, _transformation_matrix_to_super,
                              _sites, _basis)) {}

/// \brief Private constructor, implemented to initialize const members
DoFSpace::DoFSpace(
    DoFKey const &_dof_key,
    std::shared_ptr<xtal::BasicStructure const> const &_prim,
    std::optional<Eigen::Matrix3l> &&_transformation_matrix_to_super,
    std::optional<std::set<Index>> &&_sites,
    std::optional<Eigen::MatrixXd> &&_basis,
    AnisoValTraits const &_aniso_val_traits, DoFSpaceAxisInfo &&_axis_info)
    : dof_key(std::move(_dof_key)),
      is_global(_aniso_val_traits.global()),
      prim(std::move(_prim)),
      transformation_matrix_to_super(
          std::move(_transformation_matrix_to_super)),
      sites(std::move(_sites)),
      dim(get_dof_space_dimension(dof_key, *prim,
                                  transformation_matrix_to_super, sites)),
      basis(_basis.has_value() ? _basis.value()
                               : Eigen::MatrixXd::Identity(dim, dim)),
      basis_inv(basis.transpose()
                    .colPivHouseholderQr()
                    .solve(Eigen::MatrixXd::Identity(dim, dim))
                    .transpose()),
      subspace_dim(basis.cols()),
      axis_info(std::move(_axis_info)) {
  if (is_global) {
    if (transformation_matrix_to_super.has_value()) {
      std::stringstream msg;
      msg << "Error constructing DoFSpace: transformation_matrix_to_super has "
             "value for Global DoF '"
          << dof_key << "'." << std::endl;
      throw std::runtime_error(msg.str());
    }
    if (sites.has_value()) {
      std::stringstream msg;
      msg << "Error constructing DoFSpace: sites has value for Global DoF '"
          << dof_key << "'." << std::endl;
      throw std::runtime_error(msg.str());
    }
  } else {
    if (!transformation_matrix_to_super.has_value()) {
      std::stringstream msg;
      msg << "Error constructing DoFSpace: Local DoF '" << dof_key
          << "' requires transformation_matrix_to_super." << std::endl;
      throw std::runtime_error(msg.str());
    }
  }

  if (basis.rows() != dim) {
    std::stringstream msg;
    msg << "Error constructing DoFSpace: # basis rows (" << basis.rows()
        << ") != expected dimension (" << dim << ").";
    throw std::runtime_error(msg.str());
  }
  if (basis.cols() > basis.rows()) {
    std::stringstream msg;
    msg << "Error constructing DoFSpace: # basis columns (" << basis.cols()
        << ") > expected dimension (" << basis.rows() << ").";
    throw std::runtime_error(msg.str());
  }
}

/// \brief Make a DoFSpace (global or local DoF)
///
/// \param dof_key (DoFKey): A string indicating which DoF type (e.g., "disp",
///   "Hstrain", "occ")
/// \param prim (std::shared_ptr<xtal::BasicStructure const>): The prim
/// structure \param transformation_matrix_to_super
/// (std::optional<Eigen::Matrix3l>):
///   Specifies the supercell for a local DoF space. Ignored for global DoF.
/// \param sites (std::optional<std::set<Index>>): The sites included in a local
///   DoF space. Ignored for global DoF. If dof_key specifies a local DoF and
///   this does not have a value, all sites in the supercell are incldued.
/// \param basis (Eigen::MatrixXd): Allows specifying a subspace of the space
///   determined from dof_key, and for local DoF, transformation_matrix_to_super
///   and sites. The rows of `basis` correspond to prim DoF basis axes, see
///   class documentation for relations apply. If this does not have a value,
///   the standard basis (identify matrix of appropriate dimension) is used.
DoFSpace make_dof_space(
    DoFKey const &dof_key,
    std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::optional<Eigen::Matrix3l> transformation_matrix_to_super,
    std::optional<std::set<Index>> sites,
    std::optional<Eigen::MatrixXd> basis) {
  AnisoValTraits aniso_val_traits{dof_key};

  if (aniso_val_traits.global()) {
    transformation_matrix_to_super = std::nullopt;
    sites = std::nullopt;
  } else {
    if (!transformation_matrix_to_super.has_value()) {
      std::stringstream msg;
      msg << "Error in make_dof_space: Local DoF '" << dof_key
          << "' requires transformation_matrix_to_super." << std::endl;
      throw std::runtime_error(msg.str());
    }
    if (!sites.has_value()) {
      sites = std::set<Index>();
      xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
          transformation_matrix_to_super.value(), prim->basis().size());
      for (Index i = 0; i < unitcellcoord_index_converter.total_sites(); ++i) {
        sites->insert(i);
      }
    }
  }

  return DoFSpace(
      dof_key, prim, std::move(transformation_matrix_to_super),
      std::move(sites), std::move(basis), aniso_val_traits,
      DoFSpaceAxisInfo(dof_key, *prim, transformation_matrix_to_super, sites));
}

/// \brief Make a DoFSpace, convenience overload for global DoF
///
/// This overload includes only parameters necessary for a DoFSpace of global
/// DoF
DoFSpace make_global_dof_space(
    DoFKey dof_key, std::shared_ptr<xtal::BasicStructure const> const &prim,
    std::optional<Eigen::MatrixXd> basis) {
  return make_dof_space(dof_key, prim, std::nullopt, std::nullopt, basis);
}

/// \brief Make a DoFSpace, convenience overload for local DoF
///
/// This overload includes gets the `prim` and `transformation_matrix_to_super`
/// parameters necessary for a DoFSpace of global DoF from a Supercell. It will
/// also work for global DoF, though the transformation_matrix_to_super and
/// sites are not necessary for global DoF and will be ignored.
DoFSpace make_local_dof_space(
    DoFKey dof_key, std::shared_ptr<xtal::BasicStructure const> const &prim,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    std::optional<std::set<Index>> sites,
    std::optional<Eigen::MatrixXd> basis) {
  return make_dof_space(dof_key, prim, transformation_matrix_to_super, sites,
                        basis);
}

/// \brief Set DoF values from a coordinate in the DoFSpace basis
///     (continuous DoF only)
///
/// \param dof_values DoF values being set
/// \param transformation_matrix_to_super Specifies the supercell
///     the dof_values are associated with. Has no effect for global
///     DoF.
/// \param dof_space DoFSpace defining the coordinate
/// \param dof_space_coordinate The coordinate in the DoFSpace basis
///
/// \throws If `dof_space` is an occupation DoFSpace
void set_dof_value(clexulator::ConfigDoFValues &dof_values,
                   Eigen::Matrix3l const &transformation_matrix_to_super,
                   DoFSpace const &dof_space,
                   Eigen::VectorXd const &dof_space_coordinate) {
  using namespace DoFSpace_impl;
  throw_if_invalid_dof_space(transformation_matrix_to_super, dof_space);

  if (dof_space_coordinate.size() != dof_space.subspace_dim) {
    std::stringstream msg;
    msg << "Error in set_dof_value: DoF space coordinate size ("
        << dof_space_coordinate.size() << ") != # basis axes ("
        << dof_space.subspace_dim << ").";
    throw std::runtime_error(msg.str());
  }

  auto const &dof_key = dof_space.dof_key;
  auto const &basis = dof_space.basis;

  if (dof_space.is_global) {
    dof_values.global_dof_values.at(dof_key) = basis * dof_space_coordinate;
  } else {
    if (dof_key == "occ") {
      std::stringstream msg;
      msg << "Error: set_dof_value is not supported for occupation."
          << std::endl;
      throw std::runtime_error(msg.str());
    }

    auto &matrix_values = dof_values.local_dof_values.at(dof_key);
    Eigen::VectorXd vector_values = basis * dof_space_coordinate;

    auto const &axis_dof_component = dof_space.axis_info.dof_component.value();
    auto const &axis_site_index = dof_space.axis_info.site_index.value();
    for (Index i = 0; i < dof_space.dim; ++i) {
      matrix_values(axis_dof_component[i], axis_site_index[i]) =
          vector_values[i];
    }
  }
}

/// \brief Return dimension of DoFSpace
Index get_dof_space_dimension(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  if (AnisoValTraits(dof_key).global()) {
    return prim.global_dof(dof_key).dim();
  } else {
    using namespace DoFSpace_impl;
    throw_if_missing_local_dof_requirements(
        dof_key, transformation_matrix_to_super, sites);

    xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
        *transformation_matrix_to_super, prim.basis().size());
    Index dof_space_dimension = 0;
    for (Index site_index : *sites) {
      Index sublattice_index =
          unitcellcoord_index_converter(site_index).sublattice();
      xtal::Site const &site = prim.basis()[sublattice_index];
      if (dof_key == "occ") {
        dof_space_dimension += site.occupant_dof().size();
      } else if (site.has_dof(dof_key)) {
        dof_space_dimension += site.dof(dof_key).dim();
      }
    }
    return dof_space_dimension;
  }
}

/// True, if local DoF space with all sites in the supercell included
bool includes_all_sites(DoFSpace const &dof_space) {
  return dof_space.transformation_matrix_to_super.has_value() &&
         dof_space.sites.has_value() &&
         (dof_space.sites->size() ==
          dof_space.transformation_matrix_to_super->determinant() *
              dof_space.prim->basis().size());
}

/// Return true if `dof_space` is valid for `transformation_matrix_to_super`
///
/// Checks that:
/// - For local DoF, that the transformation_matrix_to_super are equivalent
bool is_valid_dof_space(Eigen::Matrix3l const &transformation_matrix_to_super,
                        DoFSpace const &dof_space) {
  if (!dof_space.is_global) {
    if (transformation_matrix_to_super !=
        dof_space.transformation_matrix_to_super.value()) {
      return false;
    }
  }
  return true;
}

/// Throw if `!is_valid_dof_space(transformation_matrix_to_super, dof_space)`
void throw_if_invalid_dof_space(
    Eigen::Matrix3l const &transformation_matrix_to_super,
    DoFSpace const &dof_space) {
  if (!is_valid_dof_space(transformation_matrix_to_super, dof_space)) {
    std::stringstream msg;
    msg << "Error: DoFSpace is not valid for given supercell." << std::endl;
    throw std::runtime_error(msg.str());
  }
}

/// \brief Constructor
///
/// \param supercell_index_converter UnitCellCoordIndexConverter for the
///     supercell in which the order parameter is begin calculated, which
///     may be different and not directly commensurate with the DoFSpace.
/// \param dof_space DoFSpace to get indexes for
///
DoFSpaceIndexConverter::DoFSpaceIndexConverter(
    xtal::UnitCellCoordIndexConverter const &_supercell_index_converter,
    DoFSpace const &dof_space)
    : prim(dof_space.prim),
      supercell_index_converter(_supercell_index_converter),
      dof_space_index_converter(
          dof_space.transformation_matrix_to_super.value(),
          dof_space.prim->basis().size()) {}

/// Perform conversion from Coordinate to DoFSpace site index
Index DoFSpaceIndexConverter::dof_space_site_index(
    xtal::Coordinate const &coord, double tol) const {
  xtal::UnitCellCoord bijk =
      xtal::UnitCellCoord::from_coordinate(*prim, coord, tol);
  return dof_space_index_converter(bijk);
}

/// Perform conversion from DoFSpace site index to supercell site index
Index DoFSpaceIndexConverter::dof_space_site_index(
    Index supercell_site_index) const {
  xtal::UnitCellCoord bijk = supercell_index_converter(supercell_site_index);
  return dof_space_index_converter(bijk);
}

/// Perform conversion from Coordinate to supercell site index
Index DoFSpaceIndexConverter::supercell_site_index(
    xtal::Coordinate const &coord, double tol) const {
  xtal::UnitCellCoord bijk =
      xtal::UnitCellCoord::from_coordinate(*prim, coord, tol);
  return supercell_index_converter(bijk);
}

/// Perform conversion from DoFSpace site index to supercell site index
Index DoFSpaceIndexConverter::supercell_site_index(
    Index dof_space_site_index) const {
  xtal::UnitCellCoord bijk = dof_space_index_converter(dof_space_site_index);
  return supercell_index_converter(bijk);
}

/// \brief Perform conversion from DoFSpace site index to supercell site index,
/// with additional translation within supercell
///
/// Equivalent to:
/// \code
/// xtal::UnitCellCoord bijk =
///     dof_space_index_converter(dof_space_site_index);
/// bijk += translation;
/// return supercell_index_converter(bijk);
/// \endcode
Index DoFSpaceIndexConverter::supercell_site_index(
    Index dof_space_site_index, xtal::UnitCell const &translation) const {
  xtal::UnitCellCoord bijk = dof_space_index_converter(dof_space_site_index);
  bijk += translation;
  return supercell_index_converter(bijk);
}

/// \brief Removes the homogeneous mode space from the DoFSpace basis
///
/// \param dof_space DoF space to remove the homogeneous mode space from.
/// Must be a DoF space for a local continuous DoF and include all sites in the
/// supercell (`includes_all_sites(dof_space)==true`), else will throw.
///
/// \returns A copy of dof_space with basis modified to remove homogeneous
/// modes.
DoFSpace exclude_homogeneous_mode_space(DoFSpace const &dof_space) {
  if (dof_space.is_global || dof_space.dof_key == "occ" ||
      !includes_all_sites(dof_space)) {
    std::stringstream msg;
    msg << "Error in exclude_homogeneous_mode_space: Must be a DoF space for a "
           "local continuous degrees of freedom that includes all sites in the "
           "supercell.";
    throw std::runtime_error(msg.str());
  }

  Eigen::MatrixXd homogeneous_mode_space =
      make_homogeneous_mode_space(dof_space);

  Eigen::MatrixXd null_space;
  if (homogeneous_mode_space.cols() == dof_space.subspace_dim) {
    null_space.resize(dof_space.dim, 0);
  } else {
    null_space = homogeneous_mode_space.transpose().fullPivLu().kernel();
  }

  return DoFSpace{dof_space.dof_key, dof_space.prim,
                  dof_space.transformation_matrix_to_super, dof_space.sites,
                  null_space};
}

/// \brief Removes the homogeneous mode space from the DoFSpace basis
///
/// \param dof_space DoF space to find the homogeneous mode space of. Should be
/// a DoF space for a local continuous DoF and include all sites in the
/// supercell.
///
/// \returns The column vector space of allowed homogeneous modes (i.e. allowed
/// rigid translations)
///
/// The prim DoF basis may not be either equal to the standard DoF basis, or
/// the same on all sublattices. The homoegeneous mode space then is limited
/// to the common part of `dof_info.basis()` for each site in the prim.
///
/// Ex: most typical, all displacements allowed on all sites:
///     prim_dof_info[0].basis(): [[dx], [dy], [dz]],
///     prim_dof_info[1].basis(): [[dx], [dy], [dz]]
///     Common basis: Rigid translations in [[dx], [dy], [dz]]
///
/// Ex: 2d displacements:
///     prim_dof_info[0].basis(): [[dx], [dy]],
///     prim_dof_info[1].basis(): [[dx], [dy]]
///     Common basis: Rigid translations are only allowed in [[dx], [dy]]
///
/// Ex: 2d displacements of differing orientation, 1d common basis:
///     prim_dof_info[0].basis(): [[dx, dy], [dz]],
///     prim_dof_info[1].basis(): [[-dx, dy], [dz]]
///     Common basis: Rigid translations are only allowed in [[dz]]
///
/// Ex: some fixed sites:
///     prim_dof_info[0].basis(): <not allowed>,
///     prim_dof_info[1].basis(): [[dx], [dy], [dz]]
///     Common basis: No rigid translations allowed
///
/// The common basis, in the standard DoF basis, is:
///     common_standard_basis = nullspace( Prod_i P(i) - I ),
/// where projector P(i), is:
///     P(i) = prim_dof_info[i].basis() * prim_dof_info[i].inv_basis()
///
/// If the common_standard_basis is not empty, the homogeneous mode space is
/// the column vector space constructed by transforming the
/// common_standard_basis back into the basis for each site in the DoFSpace:
///
///     [[ sites_dof_info[0].inv_basis() * common_standard_basis ],
///      [ sites_dof_info[1].inv_basis() * common_standard_basis ],
///      ...,
///      [ sites_dof_info[n_sites-1].inv_basis() * common_standard_basis ]]
///
/// Note that the lines above represent blocks equal to the dimension of the
/// DoF basis on each site.
Eigen::MatrixXd make_homogeneous_mode_space(DoFSpace const &dof_space) {
  if (dof_space.is_global || dof_space.dof_key == "occ" ||
      !includes_all_sites(dof_space)) {
    std::stringstream msg;
    msg << "Error in make_homogeneous_mode_space: Must be a DoF space for a "
           "local continuous degrees of freedom that includes all sites in the "
           "supercell.";
    throw std::runtime_error(msg.str());
  }

  auto const &dof_key = dof_space.dof_key;
  xtal::BasicStructure const &prim = *dof_space.prim;
  auto const &T = *dof_space.transformation_matrix_to_super;
  auto const &sites = *dof_space.sites;
  auto prim_local_dof_info = clexulator::make_local_dof_info(prim);

  /// DoFSetInfo for each sublattice
  std::vector<xtal::SiteDoFSet> prim_dof_info = prim_local_dof_info.at(dof_key);

  /// DoFSetInfo for each site in the DoFSpace with 'dof_key'
  std::vector<xtal::SiteDoFSet> sites_dof_info;
  // b: sublattice index
  // l: linear index in supercell
  // bijk: UnitCellCoord, integral site coordinates
  xtal::UnitCellCoordIndexConverter l_to_bijk(T, prim.basis().size());
  for (Index l : sites) {
    Index b = l_to_bijk(l).sublattice();
    xtal::Site const &site = prim.basis()[b];
    if (site.has_dof(dof_key)) {
      sites_dof_info.push_back(prim_dof_info[b]);
    }
  }

  // standard_dof_values = dof_info.basis() * prim_dof_values
  // dof_info.inv_basis() * standard_dof_values = prim_dof_values

  // find the common basis, in the standard DoF basis, among all sites:
  int standard_basis_dim = prim_dof_info[0].basis().rows();
  Eigen::MatrixXd I =
      Eigen::MatrixXd::Identity(standard_basis_dim, standard_basis_dim);
  Eigen::MatrixXd prod = I;
  for (auto const &sublat_dof : prim_dof_info) {
    prod = sublat_dof.basis() * sublat_dof.inv_basis() * prod;
  }
  // common_standard_basis is nullspace of (prod - I):
  Eigen::MatrixXd common_standard_basis =
      (prod - I).transpose().fullPivLu().kernel();

  if (common_standard_basis.isZero(TOL)) {
    return Eigen::MatrixXd::Zero(dof_space.dim, 0);
  }

  // construct homogeneous_mode_space by transforming common_standard_basis into
  // the site basis values, and combining for each site
  Eigen::MatrixXd homogeneous_mode_space{dof_space.dim,
                                         common_standard_basis.cols()};

  Index row = 0;  // block starting row
  Index col = 0;  // block starting column
  for (auto const &site_dof : sites_dof_info) {
    auto const &values = site_dof.inv_basis() * common_standard_basis;
    int n_rows = values.rows();
    int n_cols = values.cols();
    homogeneous_mode_space.block(row, col, n_rows, n_cols) = values;
    row += site_dof.dim();
  }

  return homogeneous_mode_space;
}

VectorSpaceMixingInfo::VectorSpaceMixingInfo(
    Eigen::MatrixXd const &column_vector_space, Eigen::MatrixXd const &subspace,
    double tol) {
  // Make a projection operator out of homogeneous mode space and project each
  // of the basis vectors onto it If they have a partial projection (not full or
  // zero) => translational modes are mixed between irreps
  Eigen::MatrixXd proj_operator = subspace * subspace.transpose();
  for (Index i = 0; i < column_vector_space.cols(); ++i) {
    Eigen::VectorXd col_projection = proj_operator * column_vector_space.col(i);
    if (col_projection.isZero(tol)) {
      axes_not_in_subspace.push_back(i);
    } else if (CASM::almost_equal(col_projection.normalized(),
                                  column_vector_space.col(i).normalized(),
                                  tol)) {
      axes_in_subspace.push_back(i);
    }

    else {
      axes_mixed_with_subspace.push_back(i);
    }
  }

  if (axes_mixed_with_subspace.size() == 0) {
    are_axes_mixed_with_subspace = false;
  }
}

}  // namespace config
}  // namespace CASM
