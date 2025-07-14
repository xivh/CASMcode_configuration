cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "Minimum OS X deployment version")

project(CASMcode_configuration VERSION 2.2.0 LANGUAGES CXX)

# set CMAKE_INSTALL_X variables
include(GNUInstallDirs)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# try to use ccache
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

##############################################
## Find dependencies

# Should find ZLIB::ZLIB
find_package(ZLIB)

# Find CASM
if(NOT DEFINED CASM_PREFIX)
  message(STATUS "CASM_PREFIX not defined")
  # try to find Python
  find_package (Python COMPONENTS Interpreter Development)
  if(DEFINED Python_EXECUTABLE)
    # if Python found, obtain CASM_PREFIX from the libcasm.casmglobal
    message(STATUS "found Python_EXECUTABLE: ${Python_EXECUTABLE}")
    message(STATUS "checking for libcasm-global")
    execute_process(
      COMMAND pip show libcasm-global
      RESULT_VARIABLE EXIT_CODE
      OUTPUT_QUIET
    )
    if (${EXIT_CODE} EQUAL 0)
      message(STATUS "found libcasm-global")
      execute_process(COMMAND ${Python_EXECUTABLE} -m libcasm.casmglobal --prefix
                      OUTPUT_VARIABLE CASM_PREFIX_RAW)
      string(STRIP ${CASM_PREFIX_RAW} CASM_PREFIX)
      message(STATUS "CASM_PREFIX: ${CASM_PREFIX}")
    else()
      message(STATUS "did not find libcasm-global")
    endif()
  endif()
endif()
if(DEFINED CASM_PREFIX)
  set(CASMcode_global_ROOT ${CASM_PREFIX}/share/CASMcode_global/cmake)
  set(CASMcode_crystallography_ROOT ${CASM_PREFIX}/share/CASMcode_crystallography/cmake)
  set(CASMcode_clexulator_ROOT ${CASM_PREFIX}/share/CASMcode_clexulator/cmake)
endif()

find_package(CASMcode_global)
if(NOT CASMcode_global_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_global")
endif()
# if successful, we have CASM::casm_global

find_package(CASMcode_crystallography)
if(NOT CASMcode_crystallography_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_crystallography")
endif()
# if successful, we have CASM::casm_crystallography

find_package(CASMcode_clexulator)
if(NOT CASMcode_clexulator_FOUND)
  message(FATAL_ERROR "CMake failed to find CASMcode_clexulator")
endif()
# if successful, we have CASM::casm_clexulator

# if no user CMAKE_INSTALL_PREFIX, use CASM_PREFIX if it exists
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED CASM_PREFIX)
    message(STATUS "CMAKE_INSTALL_PREFIX initialized to default, so updating CMAKE_INSTALL_PREFIX to CASM_PREFIX")
    set(CMAKE_INSTALL_PREFIX ${CASM_PREFIX} CACHE PATH "set CMAKE_INSTALL_PREFIX to CASM_PREFIX" FORCE)
    message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
  endif()
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

##############################################
## Build libcasm_configuration

# create libcasm_configuration
set(
  libcasm_configuration_HEADERS
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/ConfigDoFIsEquivalent.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/ConfigCompare.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/SupercellSet.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/canonical_form.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/PrimMagspinInfo.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/ConfigIsEquivalent.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/SupercellSymOp.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/PrimSymInfo.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/make_simple_structure.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/version.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/misc.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/DoFSpace_functions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/SupercellSymInfo.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/supercell_name.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/Configuration.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/config_space_analysis.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/Supercell.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/dof_space_analysis.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/copy_configuration.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/FromStructure.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/ConfigurationSet.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/Prim.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/factor_group.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/local_dof_sym_info.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/unitcellcoord_sym_info.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/occ_sym_info.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/global_dof_sym_info.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/sym_info/io/json/SymGroup_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/IrrepDecomposition.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/IrrepWedge.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/VectorSpaceSymReport.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/IrrepDecompositionImpl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/Symmetrizer.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/SimpleOrbit_impl.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/misc.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/SimpleOrbit.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/VectorSymCompare_v2.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/to_real.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/io/json/IrrepWedge_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/io/json/IrrepDecomposition_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/irreps/io/json/VectorSpaceSymReport_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccEvent.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccTrajectory.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccEventCounter.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccSystem.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/orbits.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccEventInvariants.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccEventRep.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/OccPosition.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/misc/MultiStepMethod.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/misc/LexicographicalCompare.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/io/stream/OccEvent_stream_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/io/stream/OccEventCounter_stream_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/io/json/OccEvent_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/io/json/OccSystem_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/occ_events/io/json/OccEventCounter_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/group/Group.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/group/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/group/subgroups.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/group/orbits.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/OccEventInfo.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/background_configuration.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/MakeOccEventStructures.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/perturbations.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/ConfigEnumAllOccupations.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/enumeration/ConfigurationFilter.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/io/json/Supercell_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/io/json/analysis_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/io/json/Configuration_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/ClusterSpecs.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/ClusterInvariants.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/impact_neighborhood.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/definitions.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/occ_counter.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/IntegralCluster.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/SubClusterCounter.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/orbits.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/GenericCluster.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/IntegralClusterOrbitGenerator.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/io/json/ClusterSpecs_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/io/json/EquivalentsInfo_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/io/json/IntegralClusterOrbitGenerator_json_io.hh
  ${PROJECT_SOURCE_DIR}/include/casm/configuration/clusterography/io/json/IntegralCluster_json_io.hh
)
set(
  libcasm_configuration_SOURCES
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/canonical_form.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/SupercellSet.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/PrimMagspinInfo.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/FromStructure.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/copy_configuration.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/Prim.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/ConfigurationSet.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/Supercell.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/dof_space_analysis.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/misc.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/DoFSpace_functions.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/SupercellSymInfo.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/supercell_name.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/config_space_analysis.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/Configuration.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/SupercellSymOp.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/PrimSymInfo.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/make_simple_structure.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/version.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/unitcellcoord_sym_info.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/occ_sym_info.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/local_dof_sym_info.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/factor_group.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/global_dof_sym_info.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/sym_info/io/json/SymGroup_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/VectorSpaceSymReport.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/IrrepWedge.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/IrrepDecomposition.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/VectorSymCompare_v2.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/Symmetrizer.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/IrrepDecompositionImpl.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/io/json/IrrepDecomposition_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/io/json/VectorSpaceSymReport_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/irreps/io/json/IrrepWedge_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccSystem.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccEventCounter.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccEvent.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccTrajectory.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccEventRep.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccPosition.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/OccEventInvariants.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/orbits.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/io/stream/OccEvent_stream_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/io/stream/OccEventCounter_stream_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/io/json/OccSystem_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/io/json/OccEvent_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/occ_events/io/json/OccEventCounter_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/background_configuration.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/OccEventInfo.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/ConfigurationFilter.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/perturbations.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/ConfigEnumAllOccupations.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/enumeration/MakeOccEventStructures.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/io/json/analysis_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/io/json/Supercell_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/io/json/Configuration_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/impact_neighborhood.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/ClusterSpecs.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/ClusterInvariants.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/orbits.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/occ_counter.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/IntegralCluster.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/io/json/EquivalentsInfo_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/io/json/IntegralClusterOrbitGenerator_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/io/json/ClusterSpecs_json_io.cc
  ${PROJECT_SOURCE_DIR}/src/casm/configuration/clusterography/io/json/IntegralCluster_json_io.cc
)
add_library(casm_configuration SHARED ${libcasm_configuration_SOURCES})
target_include_directories(casm_configuration
  PUBLIC
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/casm/external>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/casm/external/gzstream>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/casm/external>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/casm/external/gzstream>
)
target_compile_options(casm_configuration
  PUBLIC
    "-DCASM_CONFIGURATION_TXT_VERSION=\"${CMAKE_PROJECT_VERSION}\""
    -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=long
    -DGZSTREAM_NAMESPACE=gz
)
target_link_libraries(casm_configuration
  ZLIB::ZLIB
  ${CMAKE_DL_LIBS}
  CASM::casm_global
  CASM::casm_crystallography
  CASM::casm_clexulator
)
if(APPLE)
  set_target_properties(
    casm_configuration PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(
    casm_configuration PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()


##############################################
## Install libcasm_configuration

# install header files in <prefix>/libcasm/include/,
# while preserving directory structure
foreach ( filevar ${libcasm_configuration_HEADERS} )
  file(RELATIVE_PATH relfile ${PROJECT_SOURCE_DIR}/include/ ${filevar})
  get_filename_component( reldir ${relfile} DIRECTORY )
  install( FILES ${filevar} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${reldir} )
endforeach()

# install libcasm_configuration in <prefix>/libcasm/lib/
install(
  TARGETS casm_configuration
  EXPORT CASMcode_configurationTargets
  DESTINATION lib)

##############################################
## Python extensions

# The CMake package config and target files are installed under the Python
# package root. This is necessary to ensure that all the relative paths in the
# helloTargets.cmake resolve correctly. It also provides encapsulation.
#
# The actual path used must be selected so that consuming projects can locate it
# via `find_package`. To support finding CMake packages in the Python package
# prefix, using `find_package`s default search path of
# `<prefix>/<name>/share/<name>*/cmake/` is reasonable. Adding the Python
# package installation prefix to CMAKE_PREFIX_PATH in combination with this path
# will allow `find_package` to find this package and any other package installed
# via a Python package if the CMake and Python packages are named the same.
set(CASM_CMAKE_PACKAGE_INSTALL_SUBDIR "share/CASMcode_configuration/cmake")

install(
  EXPORT CASMcode_configurationTargets
  NAMESPACE CASM::
  DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  CASMcode_configurationConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/CASMcode_configurationConfig.cmake.in" CASMcode_configurationConfig.cmake
  INSTALL_DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

install(FILES "${PROJECT_BINARY_DIR}/CASMcode_configurationConfig.cmake"
              "${PROJECT_BINARY_DIR}/CASMcode_configurationConfigVersion.cmake"
        DESTINATION ${CASM_CMAKE_PACKAGE_INSTALL_SUBDIR})

# We are using the SKBUILD variable, which is defined when scikit-build is
# running the CMake build, to control building the Python wrapper. This allows
# the C++ project to be installed, standalone, when using the standard CMake
# build flow.
if(DEFINED SKBUILD)

  # call pybind11-config to obtain the root of the cmake package
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_ROOT_RAW)
  string(STRIP ${pybind11_ROOT_RAW} pybind11_ROOT)
  find_package(pybind11)

  # The extension modules must load:
  # - the casm_global library
  # - the casm_crystallography library
  # - the casm_clexulator library
  # - the casm_configuration library
  # They can be found by setting a relative rpath

  ### libcasm.clusterography._clusterography ###
  pybind11_add_module(_clusterography MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/clusterography.cpp")
  target_link_libraries(_clusterography PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _clusterography DESTINATION clusterography)
  if(APPLE)
    set_target_properties(
      _clusterography PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _clusterography PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.configuration._configuration ###
  pybind11_add_module(_configuration MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/configuration.cpp")
  target_link_libraries(_configuration PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _configuration DESTINATION configuration)
  if(APPLE)
    set_target_properties(
      _configuration PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _configuration PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.local_configuration._local_configuration ###
  pybind11_add_module(_local_configuration MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/local_configuration.cpp")
  target_link_libraries(_local_configuration PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _local_configuration DESTINATION local_configuration)
  if(APPLE)
    set_target_properties(
      _local_configuration PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _local_configuration PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.enumerate._enumerate ###
  pybind11_add_module(_enumerate MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/enumerate.cpp")
  target_link_libraries(_enumerate PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _enumerate DESTINATION enumerate)
  if(APPLE)
    set_target_properties(
      _enumerate PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _enumerate PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.irreps._irreps ###
  pybind11_add_module(_irreps MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/irreps.cpp")
  target_link_libraries(_irreps PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _irreps DESTINATION irreps)
  if(APPLE)
    set_target_properties(
      _irreps PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _irreps PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.occ_events._occ_events ###
  pybind11_add_module(_occ_events MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/occ_events.cpp")
  target_link_libraries(_occ_events PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _occ_events DESTINATION occ_events)
  if(APPLE)
    set_target_properties(
      _occ_events PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _occ_events PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

  ### libcasm.sym_info._sym_info ###
  pybind11_add_module(_sym_info MODULE
                      "${PROJECT_SOURCE_DIR}/python/src/sym_info.cpp")
  target_link_libraries(_sym_info PRIVATE
    CASM::casm_global
    CASM::casm_crystallography
    CASM::casm_clexulator
    casm_configuration
  )
  install(TARGETS _sym_info DESTINATION sym_info)
  if(APPLE)
    set_target_properties(
      _sym_info PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    set_target_properties(
      _sym_info PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
  endif()

endif()
