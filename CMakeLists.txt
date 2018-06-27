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
set( CMAKE_CXX_STANDARD_REQUIRED 11 CACHE STRING
   "Minimum C++ standard required for the build" )
set( CMAKE_CXX_STANDARD 14 CACHE STRING
   "C++ standard to use for the build" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL
   "(Dis)allow using compiler extensions" )

# Specify the install locations for libraries and binaries.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib )  # Needed by ROOT_GENERATE_DICTIONARY()

# Set the warning flag(s) to use.
# set( CMAKE_CXX_FLAGS "-Wall -pedantic" )


# ------------------------------------------------
# Dependencies and sub-libraries
# ------------------------------------------------

# Add ROOT system directory and require ROOT.
find_package( ROOT 6.10.06 REQUIRED MathCore HistFactory )

# Include the cmake build for the CommonStatTools submodule.
include_directories( CommonStatTools )
add_subdirectory( CommonStatTools )

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

# Handle the build of the three auxiliary classes for ATLAS-type
# labels etc. We don't want to throw this into the TRExFitter
# library object.
set( _atlas_utils_headers AtlasLabels.h AtlasStyle.h AtlasUtils.h )
set( _atlas_utils_sources AtlasLabels.C AtlasStyle.C AtlasUtils.C )
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
  TtHFitter/Common.h
  TtHFitter/ConfigParser.h
  TtHFitter/ConfigReader.h
  TtHFitter/ConfigReaderMulti.h
  TtHFitter/CorrelationMatrix.h
  TtHFitter/FitResults.h
  TtHFitter/FittingTool.h
  TtHFitter/HistoTools.h
  TtHFitter/MultiFit.h
  TtHFitter/NormFactor.h
  TtHFitter/NuisParameter.h
  TtHFitter/OverallSyst.h
  TtHFitter/Region.h
  TtHFitter/RunAsymptoticsCLs.h
  TtHFitter/RunAsymptoticsCLs_inject.h
  TtHFitter/RunSig.h
  TtHFitter/Sample.h
  TtHFitter/SampleHist.h
  TtHFitter/ShapeFactor.h
  TtHFitter/StatusLogbook.h
  TtHFitter/Systematic.h
  TtHFitter/SystematicHist.h
  TtHFitter/TtHFit.h
  TtHFitter/TthPlot.h )

# Source files for the shared/static library.
set( lib_sources
  Root/Common.C
  Root/ConfigParser.C
  Root/ConfigReader.C
  Root/ConfigReaderMulti.C
  Root/CorrelationMatrix.C
  Root/FitResults.C
  Root/FittingTool.C
  Root/HistoTools.C
  Root/MultiFit.C
  Root/NormFactor.C
  Root/NuisParameter.C
  Root/Region.C
  Root/RunAsymptoticsCLs.C
  Root/RunAsymptoticsCLs_inject.C
  Root/RunSig.C
  Root/Sample.C
  Root/SampleHist.C
  Root/ShapeFactor.C
  Root/StatusLogbook.C
  Root/Systematic.C
  Root/SystematicHist.C
  Root/TtHFit.C
  Root/TthPlot.C )

# Build the shared library.
add_library( TRExFitter SHARED ${lib_headers} ${lib_sources} )
target_include_directories( TRExFitter
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( TRExFitter CommonSystSmoothingTool ExoStats AtlasUtils ${ROOT_LIBRARIES} )
set_property( TARGET TRExFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )

# Install the libraries.
install( TARGETS TRExFitter CommonSystSmoothingTool AtlasUtils
   EXPORT TRExFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/TtHFitter )


# ------------------------------------------------
# TRExFitter executables
# ------------------------------------------------

# Helper macro for building the project's executables.
macro( TRExFitter_add_executable name )
  add_executable( ${name} ${ARGN} )
  target_link_libraries( ${name} TRExFitter CommonSystSmoothingTool ExoStats AtlasUtils )
  install( TARGETS ${name}
    EXPORT TRExFitterTargets
    RUNTIME DESTINATION bin )
endmacro( TRExFitter_add_executable )

TRExFitter_add_executable( myFit.exe util/myFit.C )


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
