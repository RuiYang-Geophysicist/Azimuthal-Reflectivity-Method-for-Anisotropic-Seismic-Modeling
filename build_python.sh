#!/bin/bash
# =============================================================================
# Build script for AzRM Python module
#
# Prerequisites:
#   - CMake 3.14+
#   - C++17 compiler (GCC 7+, Clang 5+, or MSVC 2017+)
#   - Eigen3: brew install eigen (macOS) / apt install libeigen3-dev (Ubuntu)
#   - pybind11: pip install pybind11
#   - Python 3.8+
#
# Usage:
#   ./build_python.sh          # Build release version
#   ./build_python.sh debug    # Build debug version
#   ./build_python.sh clean    # Clean build directory
#
# Author: Rui Yang, 2024
# =============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CPP_DIR="$SCRIPT_DIR/cpp"
BUILD_DIR="$CPP_DIR/build"
WEBAPP_DIR="$SCRIPT_DIR/webapp"

# Build type
BUILD_TYPE="${1:-Release}"

# Clean build
if [ "$1" == "clean" ]; then
    echo -e "${YELLOW}Cleaning build directory...${NC}"
    rm -rf "$BUILD_DIR"
    echo -e "${GREEN}Done.${NC}"
    exit 0
fi

# Debug build
if [ "$1" == "debug" ]; then
    BUILD_TYPE="Debug"
fi

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  AzRM Python Module Build${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# Check dependencies
echo -e "${YELLOW}Checking dependencies...${NC}"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}Error: Python3 not found. Please install Python 3.8+${NC}"
    exit 1
fi
PYTHON_VERSION=$(python3 --version)
echo -e "  Python: ${GREEN}$PYTHON_VERSION${NC}"

# Check CMake
if ! command -v cmake &> /dev/null; then
    echo -e "${RED}Error: CMake not found. Please install CMake 3.14+${NC}"
    exit 1
fi
CMAKE_VERSION=$(cmake --version | head -n1)
echo -e "  CMake: ${GREEN}$CMAKE_VERSION${NC}"

# Check pybind11
if ! python3 -c "import pybind11" &> /dev/null; then
    echo -e "${YELLOW}pybind11 not found. Installing...${NC}"
    pip3 install pybind11
fi
PYBIND11_VERSION=$(python3 -c "import pybind11; print(pybind11.__version__)")
echo -e "  pybind11: ${GREEN}$PYBIND11_VERSION${NC}"

# Check Eigen
EIGEN_FOUND=0
for path in /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3; do
    if [ -f "$path/Eigen/Core" ]; then
        echo -e "  Eigen3: ${GREEN}Found at $path${NC}"
        EIGEN_FOUND=1
        break
    fi
done

if [ $EIGEN_FOUND -eq 0 ]; then
    echo -e "${RED}Error: Eigen3 not found. Please install:${NC}"
    echo -e "  macOS:  brew install eigen"
    echo -e "  Ubuntu: sudo apt install libeigen3-dev"
    exit 1
fi

echo ""

# Get pybind11 CMake directory
PYBIND11_CMAKE_DIR=$(python3 -c "import pybind11; print(pybind11.get_cmake_dir())")
echo -e "${YELLOW}pybind11 CMake dir: $PYBIND11_CMAKE_DIR${NC}"

# Create build directory
echo -e "${YELLOW}Creating build directory...${NC}"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo -e "${YELLOW}Configuring with CMake...${NC}"
cmake .. \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
    -DBUILD_MEX=OFF \
    -DBUILD_PYTHON=ON \
    -DBUILD_TESTS=ON \
    -Dpybind11_DIR="$PYBIND11_CMAKE_DIR"

# Build
echo ""
echo -e "${YELLOW}Building...${NC}"
cmake --build . --config $BUILD_TYPE -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Find the built module
echo ""
echo -e "${YELLOW}Locating Python module...${NC}"

# Find the .so or .pyd file
MODULE_FILE=$(find "$BUILD_DIR" -name "azrm*.so" -o -name "azrm*.pyd" 2>/dev/null | head -n1)

if [ -z "$MODULE_FILE" ]; then
    echo -e "${RED}Error: Python module not found in build directory${NC}"
    exit 1
fi

echo -e "  Module: ${GREEN}$MODULE_FILE${NC}"

# Copy to webapp directory
echo -e "${YELLOW}Copying module to webapp directory...${NC}"
mkdir -p "$WEBAPP_DIR"
cp "$MODULE_FILE" "$WEBAPP_DIR/"

# Run tests
echo ""
echo -e "${YELLOW}Running tests...${NC}"
if [ -f "$BUILD_DIR/test_azrm" ]; then
    "$BUILD_DIR/test_azrm"
fi

# Test Python import
echo ""
echo -e "${YELLOW}Testing Python import...${NC}"
cd "$WEBAPP_DIR"
python3 -c "import azrm; print(f'azrm version: {azrm.__version__}'); print(f'OpenMP threads: {azrm.get_num_threads()}')"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Build successful!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "Python module location: ${GREEN}$WEBAPP_DIR/$(basename $MODULE_FILE)${NC}"
echo ""
echo -e "To run the webapp:"
echo -e "  cd $WEBAPP_DIR"
echo -e "  streamlit run app.py"
echo ""
