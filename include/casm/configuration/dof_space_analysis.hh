#ifndef CASM_config_dof_space_analysis
#define CASM_config_dof_space_analysis

#include "casm/clexulator/DoFSpace.hh"
#include "casm/configuration/Configuration.hh"
#include "casm/configuration/irreps/VectorSpaceSymReport.hh"

namespace CASM {
namespace clexulator {
struct DoFSpace;
}
namespace irreps {
struct VectorSpaceSymReport;
}

namespace config {

struct DoFSpaceAnalysisResults {
  DoFSpaceAnalysisResults(clexulator::DoFSpace _symmetry_adapted_dof_space,
                          irreps::VectorSpaceSymReport _symmetry_report);

  /// \brief Symmetry-adapted DoFSpace, with basis formed by
  ///     irrep decomposition
  clexulator::DoFSpace const symmetry_adapted_dof_space;

  /// \brief Summary of data associated with the action of a
  ///     symmetry group on the DoFSpace
  irreps::VectorSpaceSymReport const symmetry_report;
};

class dof_space_analysis_error : public std::runtime_error {
 public:
  dof_space_analysis_error(std::string _what) : std::runtime_error(_what) {}
  virtual ~dof_space_analysis_error() {}
};

DoFSpaceAnalysisResults dof_space_analysis(
    clexulator::DoFSpace const &dof_space, std::shared_ptr<Prim const> prim,
    std::optional<Configuration> configuration = std::nullopt,
    std::optional<bool> exclude_homogeneous_modes = std::nullopt,
    bool include_default_occ_modes = false,
    std::optional<std::map<int, int>> sublattice_index_to_default_occ =
        std::nullopt,
    std::optional<std::map<Index, int>> site_index_to_default_occ =
        std::nullopt,
    bool calc_wedges = false, std::optional<Log> log = std::nullopt);

}  // namespace config
}  // namespace CASM

#endif
