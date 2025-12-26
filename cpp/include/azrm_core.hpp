/**
 * @file azrm_core.hpp
 * @brief Azimuthal Reflectivity Method (AzRM) for Anisotropic Media
 *
 * This library implements the reflectivity method (Fryer & Frazer, 1984, 1987)
 * with numerical eigenvalue decomposition for computing reflection coefficients
 * in arbitrary anisotropic layered media.
 *
 * Supported media types:
 *   - Isotropic
 *   - VTI (Vertical Transverse Isotropy)
 *   - HTI (Horizontal Transverse Isotropy)
 *   - Orthorhombic
 *   - Monoclinic
 *   - Triclinic (fully general)
 *
 * @author Rui Yang, Tongji University, 2024
 * @author C++ implementation with OpenMP parallelization
 */

#ifndef AZRM_CORE_HPP
#define AZRM_CORE_HPP

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <vector>
#include <cmath>

namespace azrm {

// Type aliases for clarity
using Real = double;
using Complex = std::complex<double>;
using Matrix3c = Eigen::Matrix<Complex, 3, 3>;
using Matrix6c = Eigen::Matrix<Complex, 6, 6>;
using Matrix6d = Eigen::Matrix<Real, 6, 6>;
using Vector3c = Eigen::Matrix<Complex, 3, 1>;
using Vector6c = Eigen::Matrix<Complex, 6, 1>;
using Vector6d = Eigen::Matrix<Real, 6, 1>;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using MatrixXcd = Eigen::MatrixXcd;

// Constants
constexpr Real PI = 3.14159265358979323846;
constexpr Real TOL_REAL = 1e-10;
constexpr Real TOL_ABS = 1e-10;
constexpr Real REGULARIZATION_EPS = 1e-6;
constexpr Real RCOND_THRESHOLD = 1e-12;

/**
 * @struct EigenResult
 * @brief Result of eigenvalue decomposition for a single layer
 */
struct EigenResult {
    Matrix6c D;  // Eigenvector matrix (6x6)
    Matrix6c T;  // Diagonal eigenvalue matrix (6x6)
};

/**
 * @struct ReflectionCoefficients
 * @brief PP, PSv, PSh reflection coefficients for all angles
 */
struct ReflectionCoefficients {
    std::vector<Complex> Rpp;   // PP reflection coefficients
    std::vector<Complex> Rpsv;  // P-to-SV reflection coefficients
    std::vector<Complex> Rpsh;  // P-to-SH reflection coefficients
};

/**
 * @struct LayerModel
 * @brief Model parameters for layered anisotropic medium
 */
struct LayerModel {
    int nLayers;                              // Number of layers
    std::vector<Real> thickness;              // Layer thicknesses (km)
    std::vector<Real> density;                // Densities (g/cm^3)
    std::vector<Matrix6d> stiffness;          // 6x6 stiffness tensors (GPa)

    // Constructor
    LayerModel() : nLayers(0) {}
    LayerModel(int n) : nLayers(n), thickness(n), density(n), stiffness(n) {}
};

/**
 * @brief Compute eigenvalues and eigenvectors for general anisotropic media
 *
 * Uses numerical eigenvalue decomposition to solve the Christoffel equation
 * for arbitrary anisotropic media (VTI, HTI, OA, Monoclinic, Triclinic, etc.)
 *
 * @param px Horizontal slowness component in x direction (s/km)
 * @param py Horizontal slowness component in y direction (s/km)
 * @param rho Density (g/cm^3)
 * @param C 6x6 stiffness matrix in Voigt notation
 * @return EigenResult containing eigenvector matrix D and eigenvalue matrix T
 *
 * @note D columns ordered as: [up-qP, up-qS1, up-qS2, down-qP, down-qS1, down-qS2]
 *       where qP = quasi-P, qS1/qS2 = quasi-shear (sorted by |q|)
 */
EigenResult compute_eig_aniso(Real px, Real py, Real rho, const Matrix6d& C);

/**
 * @brief Compute eigenvalues and eigenvectors for VTI media (analytical)
 *
 * Uses analytical solutions by Fryer and Frazer (1987) for VTI media.
 * This is faster and more stable for T0 calculation at normal incidence.
 *
 * @param px Horizontal slowness component in x direction (s/km)
 * @param py Horizontal slowness component in y direction (s/km)
 * @param rho Density (g/cm^3)
 * @param C 6x6 stiffness matrix (VTI)
 * @return EigenResult containing eigenvector matrix D and eigenvalue matrix T
 */
EigenResult compute_eig_VTI(Real px, Real py, Real rho, const Matrix6d& C);

/**
 * @brief Compute PP, PSv, PSh reflection coefficients for anisotropic layered media
 *
 * Uses the reflectivity method (Fryer & Frazer, 1984, 1987) with numerical
 * eigenvalue decomposition. OpenMP parallelization over incidence angles.
 *
 * @param freq Frequency (Hz)
 * @param angles Incidence angles (degrees), array
 * @param model LayerModel containing layer properties
 * @param phi Azimuth angle (degrees)
 * @return ReflectionCoefficients containing Rpp, Rpsv, Rpsh at each angle
 */
ReflectionCoefficients compute_Rpp_aniso(
    Real freq,
    const std::vector<Real>& angles,
    const LayerModel& model,
    Real phi
);

/**
 * @brief Batch computation for multiple frequencies (OpenMP parallel)
 *
 * Computes reflection coefficients for multiple frequencies in parallel.
 *
 * @param frequencies Array of frequencies (Hz)
 * @param angles Incidence angles (degrees)
 * @param model LayerModel containing layer properties
 * @param phi Azimuth angle (degrees)
 * @return Vector of ReflectionCoefficients, one per frequency
 */
std::vector<ReflectionCoefficients> compute_Rpp_aniso_batch(
    const std::vector<Real>& frequencies,
    const std::vector<Real>& angles,
    const LayerModel& model,
    Real phi
);

// ============================================================================
// Utility functions
// ============================================================================

/**
 * @brief Stable matrix inversion with regularization
 */
Matrix3c stable_inv(const Matrix3c& A);

/**
 * @brief Stable linear solve A*X = B with regularization
 */
Matrix3c stable_solve(const Matrix3c& A, const Matrix3c& B);

/**
 * @brief Add regularization to near-singular matrix
 */
Matrix3c regularize_matrix(const Matrix3c& A);

/**
 * @brief Convert degrees to radians
 */
inline Real deg2rad(Real deg) { return deg * PI / 180.0; }

/**
 * @brief Convert radians to degrees
 */
inline Real rad2deg(Real rad) { return rad * 180.0 / PI; }

// ============================================================================
// C-style interface for MEX and external bindings
// ============================================================================

extern "C" {

/**
 * @brief C interface for compute_Rpp_aniso
 *
 * @param freq Frequency (Hz)
 * @param angles Array of incidence angles (degrees)
 * @param nAngles Number of angles
 * @param thickness Array of layer thicknesses (km)
 * @param density Array of densities (g/cm^3)
 * @param stiffness Flattened stiffness tensors (6*6*nLayers)
 * @param nLayers Number of layers
 * @param phi Azimuth angle (degrees)
 * @param Rpp_real Output: real part of Rpp (nAngles)
 * @param Rpp_imag Output: imaginary part of Rpp (nAngles)
 * @param Rpsv_real Output: real part of Rpsv (nAngles)
 * @param Rpsv_imag Output: imaginary part of Rpsv (nAngles)
 * @param Rpsh_real Output: real part of Rpsh (nAngles)
 * @param Rpsh_imag Output: imaginary part of Rpsh (nAngles)
 */
void azrm_compute_Rpp(
    double freq,
    const double* angles,
    int nAngles,
    const double* thickness,
    const double* density,
    const double* stiffness,  // Column-major 6x6xnLayers
    int nLayers,
    double phi,
    double* Rpp_real,
    double* Rpp_imag,
    double* Rpsv_real,
    double* Rpsv_imag,
    double* Rpsh_real,
    double* Rpsh_imag
);

/**
 * @brief C interface for batch computation over frequencies
 */
void azrm_compute_Rpp_batch(
    const double* frequencies,
    int nFreq,
    const double* angles,
    int nAngles,
    const double* thickness,
    const double* density,
    const double* stiffness,
    int nLayers,
    double phi,
    double* Rpp_real,      // nFreq x nAngles
    double* Rpp_imag,
    double* Rpsv_real,
    double* Rpsv_imag,
    double* Rpsh_real,
    double* Rpsh_imag
);

/**
 * @brief Get number of OpenMP threads
 */
int azrm_get_num_threads();

/**
 * @brief Set number of OpenMP threads
 */
void azrm_set_num_threads(int n);

} // extern "C"

} // namespace azrm

#endif // AZRM_CORE_HPP
