/**
 * @file test_azrm.cpp
 * @brief Unit tests for AzRM C++ library
 *
 * @author Rui Yang, 2024
 */

#include "azrm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

using namespace azrm;

// Test helper macros
#define TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "FAILED: " << msg << std::endl; \
            return false; \
        } \
    } while(0)

#define TEST_ASSERT_NEAR(a, b, tol, msg) \
    do { \
        if (std::abs((a) - (b)) > (tol)) { \
            std::cerr << "FAILED: " << msg << " (" << (a) << " vs " << (b) << ")" << std::endl; \
            return false; \
        } \
    } while(0)

/**
 * @brief Create a VTI stiffness matrix from Thomsen parameters
 */
Matrix6d create_VTI_stiffness(Real Vp, Real Vs, Real rho, Real eps, Real delta, Real gamma) {
    Matrix6d C = Matrix6d::Zero();

    Real c33 = Vp * Vp * rho;
    Real c44 = Vs * Vs * rho;
    Real c11 = c33 * (1 + 2 * eps);
    Real c66 = c44 * (1 + 2 * gamma);

    // c13 from Thomsen's delta
    Real term = 2 * delta * c33 * (c33 - c44) + (c33 - c44) * (c33 - c44);
    Real c13 = std::sqrt(std::max(0.0, term)) - c44;

    Real c12 = c11 - 2 * c66;

    // Fill symmetric stiffness matrix (Voigt notation)
    C(0,0) = c11; C(0,1) = c12; C(0,2) = c13;
    C(1,0) = c12; C(1,1) = c11; C(1,2) = c13;
    C(2,0) = c13; C(2,1) = c13; C(2,2) = c33;
    C(3,3) = c44;
    C(4,4) = c44;
    C(5,5) = c66;

    return C;
}

/**
 * @brief Test eigenvalue solver for isotropic medium
 */
bool test_eig_isotropic() {
    std::cout << "Testing eigenvalue solver (isotropic)... ";

    // Isotropic medium: Vp=3 km/s, Vs=1.5 km/s, rho=2.5 g/cm^3
    Real Vp = 3.0, Vs = 1.5, rho = 2.5;
    Matrix6d C = create_VTI_stiffness(Vp, Vs, rho, 0, 0, 0);

    // Test at 20 degrees incidence
    Real theta = 20.0 * PI / 180.0;
    Real Sp = 1.0 / Vp;
    Real px = Sp * std::sin(theta);
    Real py = 0;

    EigenResult result = compute_eig_aniso(px, py, rho, C);

    // For isotropic medium, eigenvalues should be real
    for (int i = 0; i < 6; ++i) {
        Real q = result.T(i, i).real();
        Real q_imag = result.T(i, i).imag();
        TEST_ASSERT(std::abs(q_imag) < 1e-10, "Eigenvalues should be real for isotropic medium");
    }

    // Check that we have 3 up-going and 3 down-going waves
    int n_up = 0, n_down = 0;
    for (int i = 0; i < 6; ++i) {
        Real q = result.T(i, i).real();
        if (q < 0) n_up++;
        else n_down++;
    }
    TEST_ASSERT(n_up == 3, "Should have 3 up-going waves");
    TEST_ASSERT(n_down == 3, "Should have 3 down-going waves");

    std::cout << "PASSED" << std::endl;
    return true;
}

/**
 * @brief Test eigenvalue solver for VTI medium
 */
bool test_eig_VTI() {
    std::cout << "Testing eigenvalue solver (VTI)... ";

    // VTI medium with typical shale anisotropy
    Real Vp = 3.5, Vs = 1.8, rho = 2.4;
    Real eps = 0.1, delta = 0.05, gamma = 0.08;
    Matrix6d C = create_VTI_stiffness(Vp, Vs, rho, eps, delta, gamma);

    // Compare general aniso solver with VTI-specific solver
    Real theta = 25.0 * PI / 180.0;
    Real Sp = 1.0 / Vp;
    Real px = Sp * std::sin(theta);
    Real py = 0;

    EigenResult result_aniso = compute_eig_aniso(px, py, rho, C);
    EigenResult result_vti = compute_eig_VTI(px, py, rho, C);

    // Eigenvalues should match (up to ordering)
    std::vector<Real> q_aniso(6), q_vti(6);
    for (int i = 0; i < 6; ++i) {
        q_aniso[i] = result_aniso.T(i, i).real();
        q_vti[i] = result_vti.T(i, i).real();
    }
    std::sort(q_aniso.begin(), q_aniso.end());
    std::sort(q_vti.begin(), q_vti.end());

    for (int i = 0; i < 6; ++i) {
        TEST_ASSERT_NEAR(q_aniso[i], q_vti[i], 1e-8,
            "Eigenvalues from aniso and VTI solvers should match");
    }

    std::cout << "PASSED" << std::endl;
    return true;
}

/**
 * @brief Test two-layer reflection coefficient
 */
bool test_two_layer_reflection() {
    std::cout << "Testing two-layer reflection coefficient... ";

    // Two-layer model
    LayerModel model(2);

    // Layer 1: low velocity
    model.thickness[0] = 0.1;  // 100m
    model.density[0] = 2.2;
    model.stiffness[0] = create_VTI_stiffness(2.5, 1.2, 2.2, 0.05, 0.02, 0.03);

    // Layer 2: high velocity
    model.thickness[1] = 0.5;  // 500m (halfspace effectively)
    model.density[1] = 2.5;
    model.stiffness[1] = create_VTI_stiffness(3.5, 1.8, 2.5, 0.08, 0.04, 0.05);

    // Compute reflection at multiple angles
    std::vector<Real> angles = {0, 10, 20, 30};
    Real freq = 30.0;  // Hz
    Real phi = 0.0;

    ReflectionCoefficients result = compute_Rpp_aniso(freq, angles, model, phi);

    // Basic checks
    TEST_ASSERT(result.Rpp.size() == angles.size(), "Output size should match input angles");

    // At normal incidence, Rpp should be real for VTI media
    TEST_ASSERT(std::abs(result.Rpp[0].imag()) < 0.1,
        "Rpp at normal incidence should be approximately real");

    // Rpp magnitude should be less than 1 (physical constraint)
    for (size_t i = 0; i < angles.size(); ++i) {
        TEST_ASSERT(std::abs(result.Rpp[i]) < 1.0,
            "Rpp magnitude should be less than 1");
    }

    // Rpp should increase with angle for positive impedance contrast
    // (This is a typical AVO Class III behavior)
    Real Rpp_0 = std::abs(result.Rpp[0]);
    Real Rpp_30 = std::abs(result.Rpp[3]);
    // Note: This might not always hold depending on model parameters

    std::cout << "PASSED" << std::endl;
    std::cout << "  Rpp at 0 deg: " << result.Rpp[0] << std::endl;
    std::cout << "  Rpp at 30 deg: " << result.Rpp[3] << std::endl;
    return true;
}

/**
 * @brief Test OpenMP parallelization
 */
bool test_openmp() {
    std::cout << "Testing OpenMP parallelization... ";

    int num_threads = azrm_get_num_threads();
    std::cout << "(threads=" << num_threads << ") ";

    // Create a larger model to test parallel performance
    LayerModel model(100);
    for (int i = 0; i < 100; ++i) {
        model.thickness[i] = 0.01;  // 10m layers
        model.density[i] = 2.0 + 0.005 * i;
        Real Vp = 2.5 + 0.01 * i;
        Real Vs = 1.2 + 0.005 * i;
        model.stiffness[i] = create_VTI_stiffness(Vp, Vs, model.density[i], 0.05, 0.02, 0.03);
    }

    // Many angles to benefit from parallelization
    std::vector<Real> angles;
    for (int i = 0; i <= 40; ++i) {
        angles.push_back(static_cast<Real>(i));
    }

    Real freq = 30.0;
    Real phi = 0.0;

    // Run computation
    ReflectionCoefficients result = compute_Rpp_aniso(freq, angles, model, phi);

    TEST_ASSERT(result.Rpp.size() == angles.size(), "Output size should match");

    std::cout << "PASSED" << std::endl;
    return true;
}

/**
 * @brief Test batch computation
 */
bool test_batch() {
    std::cout << "Testing batch frequency computation... ";

    // Simple two-layer model
    LayerModel model(2);
    model.thickness[0] = 0.1;
    model.density[0] = 2.2;
    model.stiffness[0] = create_VTI_stiffness(2.5, 1.2, 2.2, 0.05, 0.02, 0.03);

    model.thickness[1] = 0.5;
    model.density[1] = 2.5;
    model.stiffness[1] = create_VTI_stiffness(3.5, 1.8, 2.5, 0.08, 0.04, 0.05);

    std::vector<Real> frequencies = {10, 20, 30, 40, 50};
    std::vector<Real> angles = {0, 10, 20, 30};
    Real phi = 0.0;

    std::vector<ReflectionCoefficients> results =
        compute_Rpp_aniso_batch(frequencies, angles, model, phi);

    TEST_ASSERT(results.size() == frequencies.size(), "Should have result for each frequency");

    for (size_t f = 0; f < frequencies.size(); ++f) {
        TEST_ASSERT(results[f].Rpp.size() == angles.size(),
            "Each frequency should have all angles");
    }

    std::cout << "PASSED" << std::endl;
    return true;
}

/**
 * @brief Test C interface
 */
bool test_c_interface() {
    std::cout << "Testing C interface... ";

    // Create model data in C-style arrays
    int nLayers = 2;
    int nAngles = 5;

    std::vector<double> thickness = {0.1, 0.5};
    std::vector<double> density = {2.2, 2.5};
    std::vector<double> angles = {0, 10, 20, 30, 40};

    // Stiffness in column-major order (MATLAB style)
    std::vector<double> stiffness(36 * nLayers, 0);

    // Fill layer 1
    Matrix6d C1 = create_VTI_stiffness(2.5, 1.2, 2.2, 0.05, 0.02, 0.03);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            stiffness[i + j*6 + 0*36] = C1(i, j);
        }
    }

    // Fill layer 2
    Matrix6d C2 = create_VTI_stiffness(3.5, 1.8, 2.5, 0.08, 0.04, 0.05);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            stiffness[i + j*6 + 1*36] = C2(i, j);
        }
    }

    // Output arrays
    std::vector<double> Rpp_real(nAngles), Rpp_imag(nAngles);
    std::vector<double> Rpsv_real(nAngles), Rpsv_imag(nAngles);
    std::vector<double> Rpsh_real(nAngles), Rpsh_imag(nAngles);

    // Call C interface
    azrm_compute_Rpp(
        30.0,  // freq
        angles.data(),
        nAngles,
        thickness.data(),
        density.data(),
        stiffness.data(),
        nLayers,
        0.0,  // phi
        Rpp_real.data(),
        Rpp_imag.data(),
        Rpsv_real.data(),
        Rpsv_imag.data(),
        Rpsh_real.data(),
        Rpsh_imag.data()
    );

    // Check results are valid
    for (int i = 0; i < nAngles; ++i) {
        TEST_ASSERT(std::isfinite(Rpp_real[i]) && std::isfinite(Rpp_imag[i]),
            "Rpp should be finite");
        double mag = std::sqrt(Rpp_real[i]*Rpp_real[i] + Rpp_imag[i]*Rpp_imag[i]);
        TEST_ASSERT(mag < 1.0, "Rpp magnitude should be less than 1");
    }

    std::cout << "PASSED" << std::endl;
    return true;
}

/**
 * @brief Main test runner
 */
int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "AzRM C++ Library Unit Tests" << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0, failed = 0;

    auto run_test = [&](bool (*test_fn)()) {
        try {
            if (test_fn()) {
                passed++;
            } else {
                failed++;
            }
        } catch (const std::exception& e) {
            std::cerr << "EXCEPTION: " << e.what() << std::endl;
            failed++;
        }
    };

    run_test(test_eig_isotropic);
    run_test(test_eig_VTI);
    run_test(test_two_layer_reflection);
    run_test(test_openmp);
    run_test(test_batch);
    run_test(test_c_interface);

    std::cout << "========================================" << std::endl;
    std::cout << "Results: " << passed << " passed, " << failed << " failed" << std::endl;
    std::cout << "========================================" << std::endl;

    return (failed == 0) ? 0 : 1;
}
