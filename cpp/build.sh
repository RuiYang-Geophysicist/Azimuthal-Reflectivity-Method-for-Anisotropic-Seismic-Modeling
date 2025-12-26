#!/bin/bash
# ==============================================================================
# AzRM C++ Build Script
# ==============================================================================
#
# Usage:
#   ./build.sh              # Build all (MEX + Python + tests)
#   ./build.sh mex          # Build only MEX
#   ./build.sh python       # Build only Python
#   ./build.sh test         # Build and run tests
#   ./build.sh clean        # Clean build directory
#
# Prerequisites:
#   - CMake >= 3.14
#   - C++17 compiler (GCC >= 7, Clang >= 5, MSVC >= 2017)
#   - Eigen3 (brew install eigen on macOS)
#   - OpenMP (optional, for parallelization)
#   - MATLAB (optional, for MEX)
#   - Python + pybind11 (optional, pip install pybind11)
#
# ==============================================================================

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${SCRIPT_DIR}/build"
BUILD_TYPE="Release"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_header() {
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}$1${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
}

print_warning() {
    echo -e "${YELLOW}WARNING: $1${NC}"
}

print_error() {
    echo -e "${RED}ERROR: $1${NC}"
}

# Check prerequisites
check_prerequisites() {
    print_header "Checking Prerequisites"

    # CMake
    if ! command -v cmake &> /dev/null; then
        print_error "CMake not found. Please install CMake >= 3.14"
        exit 1
    fi
    echo "✓ CMake: $(cmake --version | head -n1)"

    # C++ compiler
    if command -v g++ &> /dev/null; then
        echo "✓ GCC: $(g++ --version | head -n1)"
    elif command -v clang++ &> /dev/null; then
        echo "✓ Clang: $(clang++ --version | head -n1)"
    else
        print_error "No C++ compiler found"
        exit 1
    fi

    # Eigen3
    if pkg-config --exists eigen3 2>/dev/null; then
        echo "✓ Eigen3: $(pkg-config --modversion eigen3)"
    elif [ -d "/usr/include/eigen3" ] || [ -d "/usr/local/include/eigen3" ] || [ -d "/opt/homebrew/include/eigen3" ]; then
        echo "✓ Eigen3: found in system include path"
    else
        print_warning "Eigen3 not found. Install with: brew install eigen (macOS) or apt install libeigen3-dev (Linux)"
    fi

    # OpenMP
    if [ "$(uname)" == "Darwin" ]; then
        if [ -f "/opt/homebrew/opt/libomp/lib/libomp.dylib" ] || [ -f "/usr/local/opt/libomp/lib/libomp.dylib" ]; then
            echo "✓ OpenMP: found (libomp)"
        else
            print_warning "OpenMP not found. Install with: brew install libomp"
        fi
    else
        echo "✓ OpenMP: assumed available"
    fi

    # MATLAB (optional)
    if command -v matlab &> /dev/null; then
        echo "✓ MATLAB: found"
    else
        print_warning "MATLAB not found. MEX build will be skipped."
    fi

    # Python + pybind11 (optional)
    if command -v python3 &> /dev/null; then
        echo "✓ Python: $(python3 --version)"
        if python3 -c "import pybind11" 2>/dev/null; then
            echo "✓ pybind11: found"
        else
            print_warning "pybind11 not found. Install with: pip install pybind11"
        fi
    else
        print_warning "Python3 not found. Python build will be skipped."
    fi
}

# Clean build directory
clean() {
    print_header "Cleaning Build Directory"
    rm -rf "${BUILD_DIR}"
    echo "Done."
}

# Configure and build
build() {
    local BUILD_MEX=${1:-ON}
    local BUILD_PYTHON=${2:-ON}
    local BUILD_TESTS=${3:-ON}

    print_header "Configuring Build"

    mkdir -p "${BUILD_DIR}"
    cd "${BUILD_DIR}"

    # CMake options
    CMAKE_OPTS="-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
    CMAKE_OPTS="${CMAKE_OPTS} -DBUILD_MEX=${BUILD_MEX}"
    CMAKE_OPTS="${CMAKE_OPTS} -DBUILD_PYTHON=${BUILD_PYTHON}"
    CMAKE_OPTS="${CMAKE_OPTS} -DBUILD_TESTS=${BUILD_TESTS}"
    CMAKE_OPTS="${CMAKE_OPTS} -DUSE_OPENMP=ON"

    # macOS-specific OpenMP setup
    if [ "$(uname)" == "Darwin" ]; then
        if [ -d "/opt/homebrew/opt/libomp" ]; then
            export OpenMP_ROOT="/opt/homebrew/opt/libomp"
            CMAKE_OPTS="${CMAKE_OPTS} -DOpenMP_ROOT=/opt/homebrew/opt/libomp"
        elif [ -d "/usr/local/opt/libomp" ]; then
            export OpenMP_ROOT="/usr/local/opt/libomp"
            CMAKE_OPTS="${CMAKE_OPTS} -DOpenMP_ROOT=/usr/local/opt/libomp"
        fi
    fi

    # Find MATLAB
    if [ "${BUILD_MEX}" == "ON" ]; then
        if [ -n "${MATLAB_ROOT}" ]; then
            CMAKE_OPTS="${CMAKE_OPTS} -DMatlab_ROOT_DIR=${MATLAB_ROOT}"
        elif [ -d "/Applications/MATLAB_R2024a.app" ]; then
            CMAKE_OPTS="${CMAKE_OPTS} -DMatlab_ROOT_DIR=/Applications/MATLAB_R2024a.app"
        elif [ -d "/Applications/MATLAB_R2023b.app" ]; then
            CMAKE_OPTS="${CMAKE_OPTS} -DMatlab_ROOT_DIR=/Applications/MATLAB_R2023b.app"
        fi
    fi

    echo "CMake options: ${CMAKE_OPTS}"
    cmake ${CMAKE_OPTS} ..

    print_header "Building"
    cmake --build . --parallel ${NPROC}

    echo ""
    echo -e "${GREEN}Build completed successfully!${NC}"
}

# Run tests
run_tests() {
    print_header "Running Tests"
    cd "${BUILD_DIR}"

    if [ -f "test_azrm" ]; then
        ./test_azrm
    else
        print_error "Test executable not found. Build with tests enabled first."
        exit 1
    fi
}

# Main
case "$1" in
    clean)
        clean
        ;;
    mex)
        check_prerequisites
        build ON OFF OFF
        ;;
    python)
        check_prerequisites
        build OFF ON OFF
        ;;
    test)
        check_prerequisites
        build OFF OFF ON
        run_tests
        ;;
    *)
        check_prerequisites
        build ON ON ON
        run_tests
        ;;
esac

print_header "Installation Summary"
echo "Library:     ${BUILD_DIR}/libazrm_core.a"
if [ -f "${BUILD_DIR}/azrm_mex."* ]; then
    echo "MEX file:    ${BUILD_DIR}/azrm_mex.$(ls ${BUILD_DIR}/azrm_mex.* 2>/dev/null | head -1 | sed 's/.*\.//')"
fi
if [ -f "${BUILD_DIR}/azrm."* ]; then
    echo "Python:      ${BUILD_DIR}/azrm.$(ls ${BUILD_DIR}/azrm.* 2>/dev/null | head -1 | sed 's/.*\.//')"
fi
echo ""
echo "To use in MATLAB:"
echo "  addpath('${BUILD_DIR}')"
echo "  [Rpp, Rpsv, Rpsh] = azrm_mex(freq, angles, thickness, density, stiffness, phi);"
echo ""
echo "To use in Python:"
echo "  import sys; sys.path.insert(0, '${BUILD_DIR}')"
echo "  import azrm"
echo "  Rpp, Rpsv, Rpsh = azrm.compute_Rpp(freq, angles, thickness, density, stiffness, phi)"
echo ""
