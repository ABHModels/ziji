cmake_minimum_required(VERSION 3.10)

# Set project name
project(ZijiProject VERSION 1.0 LANGUAGES CXX)

# Set the build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the required packages (GSL and OpenMP)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)

# Include directories for project-specific and external libraries (xtensor, xtl)
include_directories(
    ${PROJECT_SOURCE_DIR}/zijiray/include
    ${PROJECT_SOURCE_DIR}/external/xtensor/include
    ${PROJECT_SOURCE_DIR}/external/xtl/include
    ${GSL_INCLUDE_DIRS}    # Use GSL include directories from system
)

# Link directories for GSL (automatically handled by find_package)
# No need to specify OpenMP directories since it's handled internally by CMake

# Compiler flags (link OpenMP and GSL)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")

# Define the source files (equivalent to SRC_RAY)
file(GLOB SRC_RAY
    ${PROJECT_SOURCE_DIR}/zijiray/src/utils/*.cc
    ${PROJECT_SOURCE_DIR}/zijiray/src/GR/*.cc
    ${PROJECT_SOURCE_DIR}/zijiray/src/raytrace/*.cc
)

# Define the sources for the specific targets
set(TRANSIT_SRC ${SRC_RAY} ${PROJECT_SOURCE_DIR}/zijiray/src/transit/main.cc)
set(LPG_SRC ${SRC_RAY} ${PROJECT_SOURCE_DIR}/zijiray/src/lpg/main.cc)
set(RETURN_SRC ${SRC_RAY} ${PROJECT_SOURCE_DIR}/zijiray/src/return_rad/main.cc)

# Define the output directory for libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/ray)

# Add shared library targets (equivalent to TRANSIT_BIN, LPG_BIN, and RETURN_BIN)
add_library(transit_ray SHARED ${TRANSIT_SRC})
add_library(lpg SHARED ${LPG_SRC})
add_library(return_rad SHARED ${RETURN_SRC})

# Link the necessary libraries to the targets (GSL and OpenMP)
target_link_libraries(transit_ray PRIVATE GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX m)
target_link_libraries(lpg PRIVATE GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX m)
target_link_libraries(return_rad PRIVATE GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX m)

# Custom targets for the equivalent of lpg_bin, transit_bin, return_bin
add_custom_target(lpg_bin ALL DEPENDS lpg)
add_custom_target(transit_bin ALL DEPENDS transit_ray)
add_custom_target(return_bin ALL DEPENDS return_rad)

# Clean target to remove built binaries (equivalent to clean in Makefile)
add_custom_target(custom_clean COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_SOURCE_DIR}/bin/ray)

# Custom target to clean only CMake-generated files, not the output binaries
add_custom_target(cmake_clean
    COMMAND rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake
    COMMENT "Cleaning CMake-generated files only (CMakeCache.txt, CMakeFiles, cmake_install.cmake)"
)



# Create a custom clean target and exclude it from all builds
add_custom_target(raydata_clean
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning all files in data/data_local, except data/data_local/global_data..."
    COMMAND find data/data_local -path "data/data_local/global_data" -prune -o -type f -exec rm -f {} +
    COMMENT "Removing all ray data except global_data"
    COMMAND rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake
)

# Create a custom clean target and exclude it from all builds
add_custom_target(data_clean
    COMMAND find data -path "data/data_local/global_data" -prune -o -type f -exec rm -f {} +
    COMMENT "Removing all data except global_data"
    COMMAND rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake
)


# Ensure that all shared libraries have .so suffix, even on macOS (APPLE)
if(APPLE)
    set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")
endif()

# Force all shared libraries to use .so suffix
set_target_properties(transit_ray PROPERTIES SUFFIX ".so")
set_target_properties(lpg PROPERTIES SUFFIX ".so")
set_target_properties(return_rad PROPERTIES SUFFIX ".so")



