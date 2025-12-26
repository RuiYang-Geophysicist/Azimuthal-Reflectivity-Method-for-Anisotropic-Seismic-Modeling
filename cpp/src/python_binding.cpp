/**
 * @file python_binding.cpp
 * @brief Python bindings for AzRM using pybind11
 *
 * This module exposes the C++ AzRM library to Python using pybind11.
 * It provides a high-level Python interface for computing reflection
 * coefficients in anisotropic layered media.
 *
 * @author Rui Yang, 2024
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "azrm_core.hpp"

namespace py = pybind11;
using namespace azrm;

/**
 * @brief Compute reflection coefficients for anisotropic layered media
 *
 * @param freq Frequency (Hz)
 * @param angles NumPy array of incidence angles (degrees)
 * @param thickness NumPy array of layer thicknesses (km)
 * @param density NumPy array of densities (g/cm^3)
 * @param stiffness NumPy array of stiffness tensors, shape (nLayers, 6, 6)
 * @param phi Azimuth angle (degrees)
 * @return Tuple of (Rpp, Rpsv, Rpsh) as complex NumPy arrays
 */
py::tuple compute_reflection_coefficients(
    double freq,
    py::array_t<double> angles,
    py::array_t<double> thickness,
    py::array_t<double> density,
    py::array_t<double> stiffness,
    double phi
) {
    // Get buffer info
    auto angles_buf = angles.request();
    auto thick_buf = thickness.request();
    auto den_buf = density.request();
    auto stiff_buf = stiffness.request();

    // Validate dimensions
    if (angles_buf.ndim != 1) {
        throw std::runtime_error("angles must be 1D array");
    }
    if (thick_buf.ndim != 1) {
        throw std::runtime_error("thickness must be 1D array");
    }
    if (den_buf.ndim != 1) {
        throw std::runtime_error("density must be 1D array");
    }
    if (stiff_buf.ndim != 3) {
        throw std::runtime_error("stiffness must be 3D array with shape (nLayers, 6, 6)");
    }

    int nAngles = static_cast<int>(angles_buf.shape[0]);
    int nLayers = static_cast<int>(thick_buf.shape[0]);

    if (den_buf.shape[0] != nLayers) {
        throw std::runtime_error("density length must match number of layers");
    }
    if (stiff_buf.shape[0] != nLayers || stiff_buf.shape[1] != 6 || stiff_buf.shape[2] != 6) {
        throw std::runtime_error("stiffness must have shape (nLayers, 6, 6)");
    }

    // Build layer model
    LayerModel model(nLayers);

    double* angles_ptr = static_cast<double*>(angles_buf.ptr);
    double* thick_ptr = static_cast<double*>(thick_buf.ptr);
    double* den_ptr = static_cast<double*>(den_buf.ptr);
    double* stiff_ptr = static_cast<double*>(stiff_buf.ptr);

    std::vector<Real> angles_vec(angles_ptr, angles_ptr + nAngles);

    for (int i = 0; i < nLayers; i++) {
        model.thickness[i] = thick_ptr[i];
        model.density[i] = den_ptr[i];

        // Copy stiffness matrix (row-major in Python -> Eigen column-major)
        for (int r = 0; r < 6; r++) {
            for (int c = 0; c < 6; c++) {
                model.stiffness[i](r, c) = stiff_ptr[i * 36 + r * 6 + c];
            }
        }
    }

    // Compute reflection coefficients
    ReflectionCoefficients result = compute_Rpp_aniso(freq, angles_vec, model, phi);

    // Create output arrays
    py::array_t<std::complex<double>> Rpp(nAngles);
    py::array_t<std::complex<double>> Rpsv(nAngles);
    py::array_t<std::complex<double>> Rpsh(nAngles);

    auto Rpp_ptr = Rpp.mutable_unchecked<1>();
    auto Rpsv_ptr = Rpsv.mutable_unchecked<1>();
    auto Rpsh_ptr = Rpsh.mutable_unchecked<1>();

    for (int i = 0; i < nAngles; i++) {
        Rpp_ptr(i) = result.Rpp[i];
        Rpsv_ptr(i) = result.Rpsv[i];
        Rpsh_ptr(i) = result.Rpsh[i];
    }

    return py::make_tuple(Rpp, Rpsv, Rpsh);
}

/**
 * @brief Batch computation for multiple frequencies
 *
 * @param frequencies NumPy array of frequencies (Hz)
 * @param angles NumPy array of incidence angles (degrees)
 * @param thickness NumPy array of layer thicknesses (km)
 * @param density NumPy array of densities (g/cm^3)
 * @param stiffness NumPy array of stiffness tensors, shape (nLayers, 6, 6)
 * @param phi Azimuth angle (degrees)
 * @return Tuple of (Rpp, Rpsv, Rpsh) as complex NumPy arrays with shape (nFreq, nAngles)
 */
py::tuple compute_reflection_coefficients_batch(
    py::array_t<double> frequencies,
    py::array_t<double> angles,
    py::array_t<double> thickness,
    py::array_t<double> density,
    py::array_t<double> stiffness,
    double phi
) {
    // Get buffer info
    auto freq_buf = frequencies.request();
    auto angles_buf = angles.request();
    auto thick_buf = thickness.request();
    auto den_buf = density.request();
    auto stiff_buf = stiffness.request();

    int nFreq = static_cast<int>(freq_buf.shape[0]);
    int nAngles = static_cast<int>(angles_buf.shape[0]);
    int nLayers = static_cast<int>(thick_buf.shape[0]);

    // Build layer model
    LayerModel model(nLayers);

    double* thick_ptr = static_cast<double*>(thick_buf.ptr);
    double* den_ptr = static_cast<double*>(den_buf.ptr);
    double* stiff_ptr = static_cast<double*>(stiff_buf.ptr);

    for (int i = 0; i < nLayers; i++) {
        model.thickness[i] = thick_ptr[i];
        model.density[i] = den_ptr[i];

        for (int r = 0; r < 6; r++) {
            for (int c = 0; c < 6; c++) {
                model.stiffness[i](r, c) = stiff_ptr[i * 36 + r * 6 + c];
            }
        }
    }

    double* freq_ptr = static_cast<double*>(freq_buf.ptr);
    double* angles_ptr = static_cast<double*>(angles_buf.ptr);

    std::vector<Real> freq_vec(freq_ptr, freq_ptr + nFreq);
    std::vector<Real> angles_vec(angles_ptr, angles_ptr + nAngles);

    // Compute batch
    std::vector<ReflectionCoefficients> results =
        compute_Rpp_aniso_batch(freq_vec, angles_vec, model, phi);

    // Create output arrays
    py::array_t<std::complex<double>> Rpp({nFreq, nAngles});
    py::array_t<std::complex<double>> Rpsv({nFreq, nAngles});
    py::array_t<std::complex<double>> Rpsh({nFreq, nAngles});

    auto Rpp_ptr = Rpp.mutable_unchecked<2>();
    auto Rpsv_ptr = Rpsv.mutable_unchecked<2>();
    auto Rpsh_ptr = Rpsh.mutable_unchecked<2>();

    for (int f = 0; f < nFreq; f++) {
        for (int a = 0; a < nAngles; a++) {
            Rpp_ptr(f, a) = results[f].Rpp[a];
            Rpsv_ptr(f, a) = results[f].Rpsv[a];
            Rpsh_ptr(f, a) = results[f].Rpsh[a];
        }
    }

    return py::make_tuple(Rpp, Rpsv, Rpsh);
}

/**
 * @brief Get number of OpenMP threads
 */
int get_num_threads() {
    return azrm_get_num_threads();
}

/**
 * @brief Set number of OpenMP threads
 */
void set_num_threads(int n) {
    azrm_set_num_threads(n);
}

// =============================================================================
// Python module definition
// =============================================================================
PYBIND11_MODULE(azrm, m) {
    m.doc() = R"pbdoc(
        AzRM: Azimuthal Reflectivity Method for Anisotropic Media
        ----------------------------------------------------------

        This module provides Python bindings for the C++ AzRM library,
        which implements the reflectivity method for computing seismic
        reflection coefficients in arbitrary anisotropic layered media.

        Supported media types:
        - Isotropic
        - VTI (Vertical Transverse Isotropy)
        - HTI (Horizontal Transverse Isotropy)
        - Orthorhombic
        - Monoclinic / Triclinic

        Example usage:
            import numpy as np
            import azrm

            # Define layer model
            thickness = np.array([0.1, 0.5])  # km
            density = np.array([2.2, 2.5])    # g/cm^3
            stiffness = np.zeros((2, 6, 6))   # VTI stiffness tensors

            # Compute reflection coefficients
            angles = np.arange(0, 41)  # degrees
            Rpp, Rpsv, Rpsh = azrm.compute_Rpp(30.0, angles, thickness, density, stiffness, 0.0)
    )pbdoc";

    m.def("compute_Rpp", &compute_reflection_coefficients, R"pbdoc(
        Compute PP, PSv, PSh reflection coefficients for anisotropic layered media.

        Parameters
        ----------
        freq : float
            Frequency in Hz
        angles : ndarray
            1D array of incidence angles in degrees
        thickness : ndarray
            1D array of layer thicknesses in km
        density : ndarray
            1D array of layer densities in g/cm^3
        stiffness : ndarray
            3D array of stiffness tensors with shape (nLayers, 6, 6) in GPa
        phi : float
            Azimuth angle in degrees

        Returns
        -------
        Rpp : ndarray
            Complex PP reflection coefficients
        Rpsv : ndarray
            Complex P-to-SV reflection coefficients
        Rpsh : ndarray
            Complex P-to-SH reflection coefficients
    )pbdoc",
        py::arg("freq"),
        py::arg("angles"),
        py::arg("thickness"),
        py::arg("density"),
        py::arg("stiffness"),
        py::arg("phi") = 0.0
    );

    m.def("compute_Rpp_batch", &compute_reflection_coefficients_batch, R"pbdoc(
        Batch computation of reflection coefficients for multiple frequencies.

        Parameters
        ----------
        frequencies : ndarray
            1D array of frequencies in Hz
        angles : ndarray
            1D array of incidence angles in degrees
        thickness : ndarray
            1D array of layer thicknesses in km
        density : ndarray
            1D array of layer densities in g/cm^3
        stiffness : ndarray
            3D array of stiffness tensors with shape (nLayers, 6, 6) in GPa
        phi : float
            Azimuth angle in degrees

        Returns
        -------
        Rpp : ndarray
            Complex PP reflection coefficients with shape (nFreq, nAngles)
        Rpsv : ndarray
            Complex P-to-SV reflection coefficients
        Rpsh : ndarray
            Complex P-to-SH reflection coefficients
    )pbdoc",
        py::arg("frequencies"),
        py::arg("angles"),
        py::arg("thickness"),
        py::arg("density"),
        py::arg("stiffness"),
        py::arg("phi") = 0.0
    );

    m.def("get_num_threads", &get_num_threads, R"pbdoc(
        Get the number of OpenMP threads used for parallel computation.
    )pbdoc");

    m.def("set_num_threads", &set_num_threads, R"pbdoc(
        Set the number of OpenMP threads for parallel computation.

        Parameters
        ----------
        n : int
            Number of threads to use
    )pbdoc",
        py::arg("n")
    );

    // Version info
    m.attr("__version__") = "1.0.0";
}
