# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.1 )

# Silence some warnings on macOS with new CMake versions.
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.9 )
   cmake_policy( SET CMP0068 NEW )
endif()

# Set the project's name and version.
project( UnfoldingCode )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED 11 CACHE STRING
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
set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -O3 -g" )

# Add ROOT system directory and require ROOT.
find_package( ROOT 6.10.06 REQUIRED COMPONENTS Core MathCore Graf Hist RIO Tree Gpad )

# Public header files for the shared/static library.
set( lib_headers
    UnfoldingCode/FoldingManager.h
    UnfoldingCode/UnfoldingTools.h
    UnfoldingCode/UnfoldingResult.h )

set( lib_sources
    Root/FoldingManager.cc
    Root/UnfoldingResult.cc
    Root/UnfoldingTools.cc )

# Build the shared library.
add_library( UnfoldingCode SHARED ${lib_headers} ${lib_sources} )
target_include_directories( UnfoldingCode
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( UnfoldingCode ${ROOT_LIBRARIES} )
set_property( TARGET UnfoldingCode
   PROPERTY PUBLIC_HEADER ${lib_headers} )
target_include_directories(UnfoldingCode PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

install(TARGETS UnfoldingCode
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/UnfoldingCode )
