# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.1 )

# Silence some warnings on macOS with new CMake versions.
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.9 )
   cmake_policy( SET CMP0068 NEW )
endif()

# Set the project's name and version.
project( TRExFitter )

# Enable using CTest in the project.
include( CTest )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED 14 CACHE STRING
   "Minimum C++ standard required for the build" )
set( CMAKE_CXX_STANDARD 17 CACHE STRING
   "C++ standard to use for the build" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL
   "(Dis)allow using compiler extensions" )

# Specify the install locations for libraries and binaries.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib )  # Needed by ROOT_GENERATE_DICTIONARY()

# Set the warning flag(s) to use.
set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -O3 -g -fPIC" )


# ------------------------------------------------
# Dependencies and sub-libraries
# ------------------------------------------------

# Add ROOT system directory and require ROOT.
find_package( ROOT 6.10.06 REQUIRED COMPONENTS Core MathCore HistFactory Graf Hist RIO Tree Gpad )

# Include the cmake build for the CommonStatTools submodule.
include_directories( CommonStatTools )
if(NOT DEFINED ENV{AtlasProject})
  add_subdirectory( CommonStatTools )
endif()

set( YAML_CPP_BUILD_TESTS OFF )
set( YAML_CPP_BUILD_TOOLS OFF )

# include yaml-cpp library
include_directories( yaml-cpp )
if(NOT DEFINED ENV{AtlasProject})
  add_subdirectory( yaml-cpp )
endif()


# Handle the build of the CommonSystSmoothingTool library. That
# library comes with a cmake configuration, but introduces
# dependencies on ATLAS cmake syntax – which we are trying to
# avoid here. So ignore the shipped CMakeLists.txt, and build our
# own version of CommonSystSmoothingTool.
include_directories( CommonSystSmoothingTool )
set( _syst_smoothing_tool_sources
  CommonSystSmoothingTool/Root/PlotHist.cxx
  CommonSystSmoothingTool/Root/SmoothHist.cxx
  CommonSystSmoothingTool/Root/SmoothingTool.cxx
  CommonSystSmoothingTool/Root/TRExTools.cxx
)
add_library( CommonSystSmoothingTool SHARED ${_syst_smoothing_tool_sources} )
target_include_directories( CommonSystSmoothingTool
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( CommonSystSmoothingTool ${ROOT_LIBRARIES} )

# Add Unfolding code
include_directories( UnfoldingCode )
if(NOT DEFINED ENV{AtlasProject})
  add_subdirectory( UnfoldingCode )
endif()

# Handle the build of the three auxiliary classes for ATLAS-type
# labels etc. We don't want to throw this into the TRExFitter
# library object.
set( _atlas_utils_headers AtlasUtils/AtlasLabels.h AtlasUtils/AtlasStyle.h AtlasUtils/AtlasUtils.h )
set( _atlas_utils_sources AtlasUtils/AtlasLabels.C AtlasUtils/AtlasStyle.C AtlasUtils/AtlasUtils.C )
add_library( AtlasUtils SHARED ${_atlas_utils_headers} ${_atlas_utils_sources} )
target_include_directories( AtlasUtils
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( AtlasUtils ${ROOT_LIBRARIES} )


# ------------------------------------------------
# The actual TRExFitter library
# ------------------------------------------------

# Public header files for the shared/static library.
set( lib_headers
  TRExFitter/Common.h
  TRExFitter/ConfigParser.h
  TRExFitter/ConfigReader.h
  TRExFitter/ConfigReaderMulti.h
  TRExFitter/CorrelationMatrix.h
  TRExFitter/FitResults.h
  TRExFitter/FittingTool.h
  TRExFitter/FitUtils.h
  TRExFitter/HistoReader.h
  TRExFitter/HistoTools.h
  TRExFitter/LikelihoodScanManager.h
  TRExFitter/LimitEvent.h
  TRExFitter/MultiFit.h
  TRExFitter/NormFactor.h
  TRExFitter/NtupleReader.h
  TRExFitter/NuisParameter.h
  TRExFitter/PruningUtil.h
  TRExFitter/RankingManager.h
  TRExFitter/Region.h
  TRExFitter/Sample.h
  TRExFitter/SampleHist.h
  TRExFitter/ShapeFactor.h
  TRExFitter/StatusLogbook.h
  TRExFitter/Systematic.h
  TRExFitter/SystematicHist.h
  TRExFitter/TRExFit.h
  TRExFitter/TRExPlot.h
  TRExFitter/TruthSample.h
  TRExFitter/UnfoldingSample.h
  TRExFitter/UnfoldingSystematic.h
  TRExFitter/YamlConverter.h)

# Source files for the shared/static library.
set( lib_sources
  Root/Common.cc
  Root/ConfigParser.cc
  Root/ConfigReader.cc
  Root/ConfigReaderMulti.cc
  Root/CorrelationMatrix.cc
  Root/FitResults.cc
  Root/FittingTool.cc
  Root/FitUtils.cc
  Root/HistoReader.cc
  Root/HistoTools.cc
  Root/LikelihoodScanManager.cc
  Root/LimitEvent.cc
  Root/MultiFit.cc
  Root/NormFactor.cc
  Root/NtupleReader.cc
  Root/NuisParameter.cc
  Root/PruningUtil.cc
  Root/RankingManager.cc
  Root/Region.cc
  Root/Sample.cc
  Root/SampleHist.cc
  Root/ShapeFactor.cc
  Root/StatusLogbook.cc
  Root/Systematic.cc
  Root/SystematicHist.cc
  Root/TRExFit.cc
  Root/TRExPlot.cc
  Root/TruthSample.cc
  Root/UnfoldingSample.cc
  Root/UnfoldingSystematic.cc
  Root/YamlConverter.cc )

# Build the shared library.
add_library( TRExFitter SHARED ${lib_headers} ${lib_sources} )
target_include_directories( TRExFitter
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( TRExFitter CommonSystSmoothingTool ExoStats AtlasUtils UnfoldingCode yaml-cpp ${ROOT_LIBRARIES} )
set_property( TARGET TRExFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )
target_include_directories(TRExFitter PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

if(DEFINED ENV{AtlasProject})
  add_dependencies(TRExFitter ExoStats CommonSystSmoothingToolLib )
endif()

# Install the libraries.
install( TARGETS TRExFitter CommonSystSmoothingTool AtlasUtils
   EXPORT TRExFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/TRExFitter )


# ------------------------------------------------
# TRExFitter executables
# ------------------------------------------------

# Helper macro for building the project's executables.
macro( TRExFitter_add_executable name )
  add_executable( ${name} ${ARGN} )
  target_include_directories( ${name} PUBLIC ${ROOT_INCLUDE_DIRS} )
  target_link_libraries( ${name} TRExFitter CommonSystSmoothingTool ExoStats AtlasUtils UnfoldingCode yaml-cpp ${ROOT_LIBRARIES} )
  install( TARGETS ${name}
    EXPORT TRExFitterTargets
    RUNTIME DESTINATION bin )
endmacro( TRExFitter_add_executable )

TRExFitter_add_executable( trex-fitter util/trex-fitter.cc )

# Compile the util/hupdate.cc executable.
add_executable( hupdate.exe util/hupdate.cc )
target_include_directories( hupdate.exe PUBLIC ${ROOT_INCLUDE_DIRS} )
target_link_libraries( hupdate.exe ${ROOT_LIBRARIES} )
install( TARGETS hupdate.exe
  EXPORT TRExFitterTargets
  RUNTIME DESTINATION bin )


# ------------------------------------------------
# Install target prescription
# ------------------------------------------------

# Install the CMake description of the project. This is
# deactivated for now, because one of our dependencies (ExoStats
# aka CommonStatTools) does not define and install its targets,
# so we can't install ours either ...

# install( EXPORT TRExFitterTargets
#    FILE TRExFitterConfig-targets.cmake
#    DESTINATION cmake
#    NAMESPACE "TRExFitter::" )

# add compiler flag called TREXFITTER_HOME to not need to set environment variable any more
add_definitions( -DTREXFITTER_HOME=${CMAKE_CURRENT_SOURCE_DIR} )
