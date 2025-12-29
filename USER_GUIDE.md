# AzRM - Azimuthal Reflectivity Method for Anisotropic Media

A high-performance seismic forward modeling tool for computing synthetic seismograms in anisotropic layered media, with an interactive web-based GUI.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Algorithm](#algorithm)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [GUI User Guide](#gui-user-guide)
- [Python API](#python-api)
- [Technical Details](#technical-details)
- [References](#references)
- [Author](#author)

---

## Overview

**AzRM (Azimuthal Reflectivity Method)** implements the reflectivity method (Fryer & Frazer, 1984, 1987) with numerical eigenvalue decomposition for computing frequency-domain reflection coefficients in arbitrary anisotropic layered media. The core algorithm is written in C++ with OpenMP parallelization, exposed to Python via pybind11, and features an interactive Streamlit web GUI for easy visualization and analysis.

### Supported Media Types

| Type | Description | Required Parameters |
|------|-------------|---------------------|
| **ISO** | Isotropic | Vp, Vs, rho |
| **VTI** | Vertical Transverse Isotropy | Vp, Vs, rho, epsilon, delta, gamma |
| **HTI** | Horizontal Transverse Isotropy | Vp, Vs, rho, fracture density |
| **OA** | Orthorhombic | Vp, Vs, rho, epsilon, delta, gamma, fracture density |

---

## Features

- **High Performance**: C++ core with OpenMP parallelization over frequencies and angles
- **Full Anisotropy Support**: ISO, VTI, HTI, Orthorhombic media
- **Accurate Modeling**: Includes all internal multiples via reflectivity method
- **Interactive GUI**: Streamlit-based web interface with real-time visualization
- **Multiple Outputs**:
  - Synthetic seismic gathers (angle domain)
  - AVAZ gathers (azimuth domain)
  - AVO/AVAZ curves
  - Velocity profiles
- **Flexible Export**: PNG images (5-subplot comprehensive figures), JSON models, MAT data files (MATLAB-compatible)

---

## Algorithm

### Reflectivity Method

The algorithm solves for reflection coefficients using the propagator matrix approach:

1. **Stiffness Tensor Construction**: Build 6x6 stiffness tensors (Voigt notation) for each layer based on elastic parameters and anisotropy type

2. **Christoffel Equation**: For each layer, solve the eigenvalue problem:
   ```
   [Gamma - rho * v^2 * I] * u = 0
   ```
   where Gamma is the Christoffel matrix constructed from the stiffness tensor and slowness components

3. **Eigenvector Sorting**: Separate upgoing and downgoing waves, classify as qP (quasi-P), qS1, qS2 (quasi-shear)

4. **Recursive Reflectivity**: Apply the Kennett (1983) recursive algorithm from bottom to top:
   ```
   R_total = R_12 + T_12 * R_below * (I - R_21 * R_below)^(-1) * T_21
   ```

5. **Seismogram Synthesis**:
   - Compute Rpp(f) at all positive frequencies
   - Apply Ricker wavelet in frequency domain
   - Inverse FFT to time domain

### Key Equations

**Thomsen Parameters (VTI)**:
```
C11 = C33 * (1 + 2*epsilon)
C66 = C44 * (1 + 2*gamma)
C13 = sqrt(2*delta*C33*(C33-C44) + (C33-C44)^2) - C44
```

**HTI from Fracture Density**:
```
Using Schoenberg & Helbig (1997) penny-shaped crack model with symmetry axis rotation
```

---

## Installation

### Prerequisites

- Python 3.8+
- C++ compiler with C++17 support (clang, gcc)
- Eigen3 library
- pybind11

### Step 1: Install System Dependencies

**macOS**:
```bash
brew install eigen
```

**Windows**:
```bash
# Using vcpkg (recommended)
vcpkg install eigen3

# Or download headers manually from:
# https://eigen.tuxfamily.org/index.php?title=Main_Page
```

**Ubuntu/Debian**:
```bash
sudo apt-get install libeigen3-dev
```

### Step 2: Install Python Dependencies

```bash
cd webapp
pip install -r requirements.txt
```

The `requirements.txt` includes:
- streamlit
- numpy
- matplotlib
- pandas
- pybind11

### Step 3: Compile C++ Module

```bash
cd /path/to/00_AzRM_anisomodel_Cpp_GUI_1224
chmod +x build_python.sh
./build_python.sh
```

This will compile `azrm.cpython-*.so` in the `webapp/` directory.

### Step 4: Verify Installation

```bash
cd webapp
python3 -c "import azrm; print(f'AzRM version: {azrm.__version__}')"
```

---

## Quick Start

### Launch the GUI

```bash
cd webapp
streamlit run app.py
```

Open your browser at `http://localhost:8501`

### Minimal Python Example

```python
import numpy as np
import azrm
from stiffness import build_model_stiffness

# Define a 2-layer VTI model
layers = [
    {'type': 'VTI', 'thickness': 0.5, 'vp': 3.0, 'vs': 1.5, 'rho': 2.3,
     'epsilon': 0.05, 'delta': 0.02, 'gamma': 0.03},
    {'type': 'VTI', 'thickness': 1.0, 'vp': 4.0, 'vs': 2.0, 'rho': 2.5,
     'epsilon': 0.10, 'delta': 0.05, 'gamma': 0.08}
]

# Build stiffness tensors
thickness, density, stiffness = build_model_stiffness(layers)

# Compute Rpp at 25 Hz for angles 0-40 degrees
freq = 25.0
angles = np.arange(0, 41, 1.0)
phi = 0.0  # azimuth

Rpp, Rpsv, Rpsh = azrm.compute_Rpp(freq, angles, thickness, density, stiffness, phi)

# Note: Apply polarity fix
Rpp = -Rpp

print(f"Rpp at normal incidence: {Rpp[0]:.4f}")
```

---

## GUI User Guide

### Interface Overview

```
+------------------+---------------------------------------------+
|    SIDEBAR       |              MAIN AREA                      |
+------------------+---------------------------------------------+
| Preset Models    | Tab 1: Layer Editor                         |
| Layer Controls   |   - Edit layer properties                   |
| Forward Params   |   - Model summary table                     |
| Display Options  |                                             |
| Run Button       | Tab 2: Results                              |
|                  |   - Seismic Gather (Wiggle/Image)          |
|                  |   - AVAZ Gather (Wiggle/Image)             |
|                  |   - AVO & AVAZ Analysis                    |
|                  |   - Export Options                          |
|                  |                                             |
|                  | Tab 3: About                                |
+------------------+---------------------------------------------+
```

### Step-by-Step Tutorial

#### 1. Select or Build Your Model

**Option A: Use a Preset Model**
1. In the sidebar, find "Preset Models"
2. Select from:
   - **Two-layer VTI**: Simple VTI/VTI interface
   - **Shale-Sand**: Typical shale over sand model
   - **Fractured Reservoir**: VTI/HTI/OA model with fractures

**Option B: Build Custom Model**
1. Click **"+ Add Layer"** to add layers
2. For each layer, set:
   - **Media Type**: ISO, VTI, HTI, or OA
   - **Thickness** (km)
   - **Vp, Vs** (km/s)
   - **Density** (g/cm^3)
   - **Anisotropy parameters**:
     - VTI/OA: epsilon, delta, gamma
     - HTI/OA: fracture density

#### 2. Set Forward Modeling Parameters

In the sidebar under "Forward Modeling Parameters":

| Parameter | Description | Typical Range |
|-----------|-------------|---------------|
| theta min/max | Incidence angle range | 0-40 degrees |
| f min/max | Frequency range | 0-60 Hz |
| Angle step | Angle sampling interval | 1-2 degrees |
| Wavelet freq | Ricker wavelet dominant frequency | 20-30 Hz |
| Azimuth phi | Observation azimuth | 0-180 degrees |

#### 3. Run Forward Modeling

1. Click **"Run Forward Modeling"** button
2. Wait for computation (typically 1-3 seconds)
3. View results in the "Results" tab

#### 4. Analyze Results

**Synthetic Seismic Gather**
- Toggle between Wiggle, Image, or Both display modes
- Wiggle: Traditional variable-area display
- Image: Color-coded amplitude display

**AVAZ Seismic Gather**
- Shows amplitude variation with azimuth
- Adjust "Fixed incidence angle" slider (default: 30 degrees)
- Adjust "Azimuth step" for resolution
- Useful for detecting HTI/OA anisotropy

**AVO & AVAZ Analysis**
- Select interface for analysis
- **AVO Curve**: |Rpp| vs incidence angle
  - Adjust azimuth angle slider (0-180°) to observe azimuthal variation
  - Useful for analyzing AVO response at different azimuths
- **AVAZ Curve**: |Rpp| vs azimuth
  - Adjust fixed incidence angle slider to observe azimuthal anisotropy
  - Useful for detecting HTI/OA anisotropy patterns
- Anisotropy indicator shows variation percentage

#### 5. Export Results

Three export formats available:

- **PNG**: High-resolution comprehensive figure (300 dpi) containing 5 subplots:
  1. Velocity profile
  2. AVO wiggle display
  3. AVO image display
  4. AVAZ wiggle display
  5. AVAZ image display

- **JSON**: Model parameters for reproducibility and model sharing

- **MAT**: MATLAB-compatible data file containing:
  - AVO data: time, angles, seismogram, reflection coefficients (real/imaginary)
  - AVAZ data: time, azimuths, seismogram, incidence angle
  - Model parameters: layer properties, frequency settings
  - All data formatted for easy import into MATLAB

### Tips for Best Results

1. **For VTI models**: Ensure |epsilon| < 0.3, |delta| < 0.3, |gamma| < 0.3
2. **For HTI models**: Fracture density typically 0.01-0.1
3. **For strong AVAZ effect**: Use HTI or OA layers with phi != 0 or 90 degrees
4. **For stable computation**: Avoid very thin layers (< 0.01 km)

---

## Python API

### Core Functions

#### `azrm.compute_Rpp(freq, angles, thickness, density, stiffness, phi)`

Compute reflection coefficients at a single frequency.

**Parameters**:
- `freq` (float): Frequency in Hz
- `angles` (ndarray): Incidence angles in degrees
- `thickness` (list): Layer thicknesses in km
- `density` (list): Layer densities in g/cm^3
- `stiffness` (list of 6x6 arrays): Stiffness tensors in GPa
- `phi` (float): Azimuth angle in degrees

**Returns**:
- `Rpp` (ndarray): Complex PP reflection coefficients
- `Rpsv` (ndarray): Complex P-SV reflection coefficients
- `Rpsh` (ndarray): Complex P-SH reflection coefficients

**Note**: Apply polarity fix with `Rpp = -Rpp`

#### `azrm.compute_Rpp_batch(frequencies, angles, thickness, density, stiffness, phi)`

Batch computation for multiple frequencies (parallelized).

**Parameters**:
- `frequencies` (ndarray): Array of frequencies in Hz
- (other parameters same as `compute_Rpp`)

**Returns**:
- `Rpp` (ndarray): Shape (nFreq, nAngles), complex
- `Rpsv`, `Rpsh`: Same shape

#### `azrm.get_num_threads()` / `azrm.set_num_threads(n)`

Get/set OpenMP thread count.

### Stiffness Module

```python
from stiffness import (
    build_model_stiffness,      # Build full model
    stiffness_isotropic,        # ISO stiffness tensor
    stiffness_vti,              # VTI stiffness tensor
    stiffness_hti,              # HTI stiffness tensor
    stiffness_orthorhombic,     # OA stiffness tensor
    preset_two_layer_vti,       # Preset model
    preset_shale_sand,          # Preset model
    preset_fractured_reservoir  # Preset model
)
```

---

## Technical Details

### Project Structure

```
00_AzRM_anisomodel_Cpp_GUI_1224/
├── cpp/                          # C++ core library
│   ├── include/
│   │   └── azrm_core.hpp         # Header with API declarations
│   ├── src/
│   │   ├── azrm_core.cpp         # Core reflectivity algorithm
│   │   ├── python_binding.cpp    # pybind11 Python bindings
│   │   └── mex_interface.cpp     # MATLAB MEX interface
│   ├── test/
│   │   └── test_azrm.cpp         # Unit tests
│   ├── CMakeLists.txt
│   └── build.sh
├── webapp/                       # Streamlit GUI
│   ├── app.py                    # Main application (~900 lines)
│   ├── stiffness.py              # Stiffness tensor construction
│   ├── seismic_plot.py           # Visualization functions
│   └── requirements.txt          # Python dependencies
├── build_python.sh               # Build script for Python module
├── *.m                           # Original MATLAB implementations
└── README.md                     # This file
```

### Performance

| Model Size | Frequencies | Angles | Time (8 threads) |
|------------|-------------|--------|------------------|
| 2 layers   | 256         | 41     | ~0.3 s           |
| 5 layers   | 512         | 41     | ~1.2 s           |
| 10 layers  | 512         | 41     | ~3.5 s           |

### Numerical Stability

- Regularization for near-singular matrices (epsilon = 1e-6)
- Condition number checking (threshold = 1e-12)
- Phase continuity correction at normal incidence

---

## References

1. **Yang R, Chen H, Guo Z, et al.** An effective azimuthal reflectivity modeling (AzRM) tool for generating seismic data in anisotropic shale reservoirs. *Geophysics*, 2025, 90(5): 1-76.

2. **Fryer G J, Frazer L N.** Seismic waves in stratified anisotropic media. *Geophysical Journal International*, 1984, 78(3): 691-710.

3. **Fryer G J, Frazer L N.** Seismic waves in stratified anisotropic media—II. Elastodynamic eigensolutions for some anisotropic systems. *Geophysical Journal International*, 1987, 91(1): 73-101.

4. **Kennett B L.** Seismic wave propagation in stratified media. *Cambridge University Press*, 1983.

5. **Thomsen L.** Weak elastic anisotropy. *Geophysics*, 1986, 51(10): 1954-1966.

6. **Schoenberg M, Helbig K.** Orthorhombic media: Modeling elastic wave behavior in a vertically fractured earth. *Geophysics*, 1997, 62(6): 1954-1974.

---

## Author

**Rui Yang**
Tongji University / Stanford University
2024-2025

---

## License

This software is provided for academic and research purposes. Please cite the relevant publications when using this tool in your research.

---

## Troubleshooting

### C++ module not found

```bash
# Rebuild the module
./build_python.sh

# Check if the .so file exists
ls webapp/*.so
```

### Streamlit won't start

```bash
# Use explicit Python module invocation
python3 -m streamlit run app.py --server.headless true
```

### Numerical instability warnings

- Try reducing layer contrast (velocity jumps)
- Increase layer thickness (avoid < 0.01 km)
- Check for physically unrealistic parameters

### AVAZ shows no variation

- AVAZ variation requires HTI or OA layers
- Ensure fracture_density > 0 for HTI layers
- Try different azimuth angles
