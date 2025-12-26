cmake_minimum_required(VERSION 3.14)
project(AzRM VERSION 1.0.0 LANGUAGES CXX)

# ==============================================================================
# Build options
# ==============================================================================
option(BUILD_MEX "Build MATLAB MEX interface" ON)
option(BUILD_PYTHON "Build Python bindings" ON)
option(BUILD_TESTS "Build tests" ON)
option(USE_OPENMP "Enable OpenMP parallelization" ON)

# ==============================================================================
# C++ standard
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    # Enable architecture-specific optimizations
    add_compile_options(-march=native)
endif()

# ==============================================================================
# Dependencies
# ==============================================================================

# Eigen (header-only)
find_package(Eigen3 CONFIG)
if(NOT Eigen3_FOUND)
    # Fallback: try to find Eigen headers directly
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core
        HINTS /usr/include/eigen3 /usr/local/include/eigen3
              /opt/homebrew/include/eigen3)
    if(EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3 headers at: ${EIGEN3_INCLUDE_DIR}")
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR "Eigen3 not found. Please install Eigen3.")
    endif()
endif()

# OpenMP
if(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    else()
        message(WARNING "OpenMP not found. Building without parallelization.")
        set(USE_OPENMP OFF)
    endif()
endif()

# ==============================================================================
# Core library
# ==============================================================================
add_library(azrm_core STATIC
    src/azrm_core.cpp
)

target_include_directories(azrm_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(azrm_core PUBLIC
    Eigen3::Eigen
)

if(USE_OPENMP AND OpenMP_CXX_FOUND)
    target_link_libraries(azrm_core PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(azrm_core PUBLIC USE_OPENMP)
endif()

# ==============================================================================
# MATLAB MEX interface
# ==============================================================================
if(BUILD_MEX)
    # Find MATLAB
    find_package(Matlab COMPONENTS MX_LIBRARY MEX_COMPILER)

    if(Matlab_FOUND)
        message(STATUS "MATLAB found: ${Matlab_ROOT_DIR}")

        matlab_add_mex(
            NAME azrm_mex
            SRC src/mex_interface.cpp
            LINK_TO azrm_core
        )

        # Additional MEX compile flags
        target_compile_definitions(azrm_mex PRIVATE
            MATLAB_MEX_FILE
        )

        # Install MEX file
        install(TARGETS azrm_mex
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/matlab
        )
    else()
        message(WARNING "MATLAB not found. Skipping MEX build.")
        message(STATUS "To build MEX, set Matlab_ROOT_DIR or add MATLAB to PATH")
    endif()
endif()

# ==============================================================================
# Python bindings
# ==============================================================================
if(BUILD_PYTHON)
    # Find Python and pybind11
    find_package(Python3 COMPONENTS Interpreter Development)
    find_package(pybind11 CONFIG)

    if(Python3_FOUND AND pybind11_FOUND)
        message(STATUS "Python found: ${Python3_EXECUTABLE}")
        message(STATUS "pybind11 found: ${pybind11_VERSION}")

        pybind11_add_module(azrm MODULE
            src/python_binding.cpp
        )

        target_link_libraries(azrm PRIVATE azrm_core)

        # Install Python module
        install(TARGETS azrm
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/python
        )
    else()
        message(WARNING "Python or pybind11 not found. Skipping Python build.")
        message(STATUS "Install pybind11: pip install pybind11")
    endif()
endif()

# ==============================================================================
# Tests
# ==============================================================================
if(BUILD_TESTS)
    enable_testing()

    add_executable(test_azrm
        test/test_azrm.cpp
    )

    target_link_libraries(test_azrm PRIVATE azrm_core)

    add_test(NAME test_azrm COMMAND test_azrm)
endif()

# ==============================================================================
# Install
# ==============================================================================
install(TARGETS azrm_core
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

install(FILES include/azrm_core.hpp
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

# ==============================================================================
# Summary
# ==============================================================================
message(STATUS "")
message(STATUS "=== AzRM Build Configuration ===")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "USE_OPENMP: ${USE_OPENMP}")
message(STATUS "BUILD_MEX: ${BUILD_MEX}")
message(STATUS "BUILD_PYTHON: ${BUILD_PYTHON}")
message(STATUS "BUILD_TESTS: ${BUILD_TESTS}")
message(STATUS "================================")
message(STATUS "")
