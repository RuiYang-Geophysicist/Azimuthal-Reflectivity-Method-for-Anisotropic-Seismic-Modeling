/**
 * @file azrm_core.cpp
 * @brief Implementation of AzRM core algorithms
 *
 * @author Rui Yang, Tongji University, 2024
 * @author C++ implementation with OpenMP parallelization
 */

#include "azrm_core.hpp"
#include <algorithm>
#include <numeric>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace azrm {

// ============================================================================
// Utility functions implementation
// ============================================================================

Matrix3c regularize_matrix(const Matrix3c& A) {
    // Compute condition number estimate
    Eigen::JacobiSVD<Matrix3c> svd(A);
    auto singularValues = svd.singularValues();

    Real maxSV = singularValues(0);
    Real minSV = singularValues(2);

    if (maxSV < 1e-20) {
        // Matrix is essentially zero, add identity
        return A + Complex(REGULARIZATION_EPS, 0) * Matrix3c::Identity();
    }

    Real rcond = minSV / maxSV;

    if (!std::isfinite(rcond) || rcond < RCOND_THRESHOLD) {
        // Add regularization proportional to Frobenius norm
        Real frobNorm = A.norm();
        return A + Complex(REGULARIZATION_EPS * frobNorm, 0) * Matrix3c::Identity();
    }

    return A;
}

Matrix3c stable_inv(const Matrix3c& A) {
    // Regularize if needed
    Matrix3c Areg = regularize_matrix(A);

    // Check for invalid values
    bool hasInvalid = false;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!std::isfinite(Areg(i,j).real()) || !std::isfinite(Areg(i,j).imag())) {
                hasInvalid = true;
                Areg(i,j) = Complex(0, 0);
            }
        }
    }
    if (hasInvalid) {
        Areg += Complex(REGULARIZATION_EPS, 0) * Matrix3c::Identity();
    }

    // Try LU decomposition
    Eigen::FullPivLU<Matrix3c> lu(Areg);

    if (lu.isInvertible()) {
        Matrix3c result = lu.inverse();

        // Check result validity
        bool resultValid = true;
        for (int i = 0; i < 3 && resultValid; ++i) {
            for (int j = 0; j < 3 && resultValid; ++j) {
                if (!std::isfinite(result(i,j).real()) || !std::isfinite(result(i,j).imag())) {
                    resultValid = false;
                }
            }
        }

        if (resultValid) {
            return result;
        }
    }

    // Fallback to pseudo-inverse using SVD
    Eigen::JacobiSVD<Matrix3c> svd(Areg, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto singularValues = svd.singularValues();

    // Compute pseudo-inverse
    Matrix3c Sinv = Matrix3c::Zero();
    Real tolerance = 1e-10 * singularValues(0);
    for (int i = 0; i < 3; ++i) {
        if (singularValues(i) > tolerance) {
            Sinv(i,i) = Complex(1.0 / singularValues(i), 0);
        }
    }

    return svd.matrixV() * Sinv * svd.matrixU().adjoint();
}

Matrix3c stable_solve(const Matrix3c& A, const Matrix3c& B) {
    // Regularize if needed
    Matrix3c Areg = regularize_matrix(A);

    // Check for invalid values
    bool hasInvalid = false;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (!std::isfinite(Areg(i,j).real()) || !std::isfinite(Areg(i,j).imag())) {
                hasInvalid = true;
                Areg(i,j) = Complex(0, 0);
            }
        }
    }
    if (hasInvalid) {
        Areg += Complex(REGULARIZATION_EPS, 0) * Matrix3c::Identity();
    }

    // Solve using LU decomposition
    Eigen::FullPivLU<Matrix3c> lu(Areg);
    Matrix3c result = lu.solve(B);

    // Check result validity
    bool resultValid = true;
    for (int i = 0; i < 3 && resultValid; ++i) {
        for (int j = 0; j < 3 && resultValid; ++j) {
            if (!std::isfinite(result(i,j).real()) || !std::isfinite(result(i,j).imag())) {
                resultValid = false;
            }
        }
    }

    if (!resultValid) {
        // Fallback: X = pinv(A) * B
        Matrix3c Ainv = stable_inv(Areg);
        result = Ainv * B;
    }

    return result;
}

// ============================================================================
// compute_eig_aniso implementation
// ============================================================================

namespace {

/**
 * @brief Reorder eigenpairs for general anisotropic media
 *
 * Convention:
 *   - First 3 columns: up-going waves (Re(q) < 0 or energy flux upward)
 *   - Last 3 columns: down-going waves (Re(q) > 0 or energy flux downward)
 *   - Within each group: sorted by |q| ascending (fastest to slowest)
 */
void reorder_eigenpairs_general(
    const Eigen::Matrix<Complex, 6, 6>& Din,
    const Eigen::Matrix<Complex, 6, 6>& Tin,
    Eigen::Matrix<Complex, 6, 6>& Dout,
    Eigen::Matrix<Complex, 6, 6>& Tout
) {
    // Extract eigenvalues
    Vector6c q;
    for (int i = 0; i < 6; ++i) {
        q(i) = Tin(i, i);
    }

    // Step 1: Compute energy flux for each mode
    Eigen::Matrix<Real, 6, 1> energy_flux_z;
    for (int j = 0; j < 6; ++j) {
        Vector3c u = Din.block<3,1>(0, j);    // displacement
        Vector3c tau = Din.block<3,1>(3, j);  // traction
        // Vertical energy flux - use explicit evaluation
        Complex flux = (tau.adjoint() * u).value() + (u.adjoint() * tau).value();
        energy_flux_z(j) = flux.real();
    }

    // Step 2: Classify as up-going or down-going
    std::vector<int> idx_up, idx_down;

    for (int j = 0; j < 6; ++j) {
        if (q(j).real() > TOL_REAL) {
            idx_down.push_back(j);  // Positive Re(q) -> down-going
        } else if (q(j).real() < -TOL_REAL) {
            idx_up.push_back(j);    // Negative Re(q) -> up-going
        } else {
            // Purely imaginary or zero: use energy flux
            if (energy_flux_z(j) > 0) {
                idx_down.push_back(j);
            } else {
                idx_up.push_back(j);
            }
        }
    }

    // Fallback: if not exactly 3+3, force split by Re(q)
    if (idx_up.size() != 3 || idx_down.size() != 3) {
        std::vector<std::pair<Real, int>> sorted_q(6);
        for (int i = 0; i < 6; ++i) {
            sorted_q[i] = {q(i).real(), i};
        }
        std::sort(sorted_q.begin(), sorted_q.end());

        idx_up.clear();
        idx_down.clear();
        for (int i = 0; i < 3; ++i) {
            idx_up.push_back(sorted_q[i].second);
        }
        for (int i = 3; i < 6; ++i) {
            idx_down.push_back(sorted_q[i].second);
        }
    }

    // Step 3: Sort each group by |q| ascending
    auto sort_by_abs_q = [&q](int a, int b) {
        return std::abs(q(a)) < std::abs(q(b));
    };

    std::sort(idx_up.begin(), idx_up.end(), sort_by_abs_q);
    std::sort(idx_down.begin(), idx_down.end(), sort_by_abs_q);

    // Final ordering: [up-fast, up-mid, up-slow, down-fast, down-mid, down-slow]
    std::vector<int> final_order;
    for (int i : idx_up) final_order.push_back(i);
    for (int i : idx_down) final_order.push_back(i);

    // Build output matrices
    for (int j = 0; j < 6; ++j) {
        Dout.col(j) = Din.col(final_order[j]);
        Tout(j, j) = q(final_order[j]);
    }

    // Step 5: Normalize eigenvectors using energy flux
    for (int j = 0; j < 6; ++j) {
        Vector3c u = Dout.block<3,1>(0, j);
        Vector3c tau = Dout.block<3,1>(3, j);

        Complex energy_flux = (tau.adjoint() * u).value() + (u.adjoint() * tau).value();

        if (std::abs(energy_flux) > 1e-14) {
            // Energy normalization: |flux| = 1
            Real epsilon = 1.0 / std::sqrt(std::abs(energy_flux));
            Dout.col(j) *= epsilon;
        } else {
            // Degenerate case: use Euclidean norm
            Dout.col(j) /= Dout.col(j).norm();
        }

        // Sign convention compatible with VTI analytical code
        Vector3c u_norm = Dout.block<3,1>(0, j);
        Real ux = u_norm(0).real();
        Real uy = u_norm(1).real();
        Real uz = u_norm(2).real();

        if (std::abs(ux) > TOL_ABS) {
            if (ux < 0) Dout.col(j) *= -1.0;
        } else if (std::abs(uy) > TOL_ABS) {
            if (uy < 0) Dout.col(j) *= -1.0;
        } else {
            if (uz < 0) Dout.col(j) *= -1.0;
        }
    }
}

} // anonymous namespace

EigenResult compute_eig_aniso(Real px, Real py, Real rho, const Matrix6d& C) {
    EigenResult result;
    result.D = Matrix6c::Zero();
    result.T = Matrix6c::Zero();

    // Extract stiffness components (Triclinic - most general)
    Real c11 = C(0,0), c12 = C(0,1), c13 = C(0,2), c16 = C(0,5);
    Real c22 = C(1,1), c23 = C(1,2), c26 = C(1,5);
    Real c33 = C(2,2), c36 = C(2,5);
    Real c44 = C(3,3), c45 = C(3,4);
    Real c55 = C(4,4);
    Real c66 = C(5,5);
    // Note: c14, c15, c24, c25, c34, c35, c46, c56 are zero for VTI and many media types

    // Build the 6x6 system matrix A
    // State vector: [ux, uy, uz, tau_xz, tau_yz, tau_zz]^T

    // S matrix elements (stress-displacement relation in horizontal plane)
    Real s11 = px*px*(c11 - c13*c13/c33) + 2*px*py*(c16 - c13*c36/c33)
               + py*py*(c66 - c36*c36/c33) - rho;
    Real s12 = px*px*(c16 - c13*c36/c33)
               + px*py*(c12 + c66 - c13*c23/c33 - c36*c36/c33)
               + py*py*(c26 - c23*c36/c33);
    Real s22 = px*px*(c66 - c36*c36/c33) + 2*px*py*(c26 - c23*c36/c33)
               + py*py*(c22 - c23*c23/c33) - rho;

    // T matrix (relates displacement to traction gradient)
    Eigen::Matrix<Real, 3, 3> Ta;
    Ta << 0,                           0,                           -px,
          0,                           0,                           -py,
          -(px*c13 + py*c36)/c33,     -(px*c36 + py*c23)/c33,       0;

    // C matrix (compliance-related)
    Real det_c45 = c44*c55 - c45*c45;
    if (std::abs(det_c45) < 1e-20) {
        det_c45 = (det_c45 >= 0 ? 1.0 : -1.0) * 1e-20;  // Regularization
    }

    Eigen::Matrix<Real, 3, 3> Ca;
    Ca << -c44/det_c45,   c45/det_c45,  0,
           c45/det_c45,  -c55/det_c45,  0,
           0,             0,           -1/c33;

    // S matrix
    Eigen::Matrix<Real, 3, 3> Sa;
    Sa << s11,  s12,  0,
          s12,  s22,  0,
          0,    0,   -rho;

    // Assemble full 6x6 system matrix
    Matrix6d A;
    A.block<3,3>(0,0) = Ta;
    A.block<3,3>(0,3) = Ca;
    A.block<3,3>(3,0) = Sa;
    A.block<3,3>(3,3) = Ta.transpose();

    // Solve eigenvalue problem
    Eigen::EigenSolver<Matrix6d> eigensolver(A);

    Matrix6c D_raw = eigensolver.eigenvectors();
    Vector6c eigenvalues = eigensolver.eigenvalues();

    // Create diagonal eigenvalue matrix
    Matrix6c T_raw = Matrix6c::Zero();
    for (int i = 0; i < 6; ++i) {
        T_raw(i, i) = eigenvalues(i);
    }

    // Reorder eigenpairs
    reorder_eigenpairs_general(D_raw, T_raw, result.D, result.T);

    return result;
}

// ============================================================================
// compute_eig_VTI implementation (analytical solution)
// ============================================================================

EigenResult compute_eig_VTI(Real px, Real py, Real rho, const Matrix6d& C) {
    EigenResult result;
    result.D = Matrix6c::Zero();
    result.T = Matrix6c::Zero();

    // Extract stiffness coefficients for VTI medium
    Real c11 = C(0,0), c12 = C(0,1), c13 = C(0,2);
    Real c33 = C(2,2), c44 = C(3,3), c66 = C(5,5);

    // Auxiliary scalar terms
    Real s11 = px*px * (c11 - c13*c13/c33) + py*py * c66 - rho;
    Real s12 = px*py * (c12 + c66 - c13*c13/c33);
    Real s22 = px*px * c66 + py*py * (c11 - c13*c13/c33) - rho;

    // Construct system matrix A (Fryer & Frazer 1987, Eq. 4.9)
    Matrix6d A;
    A << 0, 0, -px*c13/c33, s11, s12, 0,
         0, 0, -py*c13/c33, s12, s22, 0,
         -px, -py, 0, 0, 0, -rho,
         -1/c44, 0, 0, 0, 0, -px,
         0, -1/c44, 0, 0, 0, -py,
         0, 0, -1/c33, -px*c13/c33, -py*c13/c33, 0;

    // Solve eigenproblem
    Eigen::EigenSolver<Matrix6d> eigensolver(A);
    Vector6c eigenvalues = eigensolver.eigenvalues();

    // Extract real parts and sort
    std::vector<std::pair<Real, int>> qi_sorted(6);
    for (int i = 0; i < 6; ++i) {
        qi_sorted[i] = {eigenvalues(i).real(), i};
    }
    std::sort(qi_sorted.begin(), qi_sorted.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Reorder: [qD_T > qD_S > qD_P > 0 > qU_P > qU_S > qU_T]
    // Output order: [up-P, up-SV, up-SH, down-P, down-SV, down-SH]
    // = [qi[3], qi[4], qi[5], qi[2], qi[1], qi[0]] in sorted order
    Real q[6];
    q[0] = qi_sorted[3].first;  // up-P
    q[1] = qi_sorted[4].first;  // up-SV
    q[2] = qi_sorted[5].first;  // up-SH
    q[3] = qi_sorted[2].first;  // down-P
    q[4] = qi_sorted[1].first;  // down-SV
    q[5] = qi_sorted[0].first;  // down-SH

    // Set eigenvalue matrix
    for (int i = 0; i < 6; ++i) {
        result.T(i, i) = Complex(q[i], 0);
    }

    // Compute eigenvectors analytically
    for (int i = 0; i < 6; ++i) {
        Vector6c eigvec;
        Real qi = q[i];

        if (i == 0 || i == 3) {
            // P-wave modes
            Real denom = (c13 + c44) * qi;
            if (std::abs(denom) < 1e-20) denom = 1e-20;
            Real L_P = -(c11 * px*px + c44 * qi*qi - rho) / denom;

            Real ux = px, uy = py, uz = L_P;
            Real tau1 = -c44 * px * (L_P + qi);
            Real tau2 = -c44 * py * (L_P + qi);
            Real tau3 = -(c13 * px*px + c33 * qi * L_P);

            Vector3c uP(ux, uy, uz);
            Vector3c tauP(tau1, tau2, tau3);

            Complex flux = (tauP.adjoint() * uP).value() + (uP.adjoint() * tauP).value();
            Real epsilon_p = 1.0 / std::sqrt(std::abs(flux) + 1e-30);

            eigvec.head<3>() = epsilon_p * uP;
            eigvec.tail<3>() = epsilon_p * tauP;

        } else if (i == 1 || i == 4) {
            // SV-wave modes
            Real denom = (c13 + c44) * qi;
            if (std::abs(denom) < 1e-20) denom = 1e-20;
            Real L_SV = -(c11 * px*px + c44 * qi*qi - rho) / denom;

            Real ux = px, uy = py, uz = L_SV;
            Real tau1 = -c44 * px * (L_SV + qi);
            Real tau2 = -c44 * py * (L_SV + qi);
            Real tau3 = -(c13 * px*px + c33 * qi * L_SV);

            Vector3c uSV(ux, uy, uz);
            Vector3c tauSV(tau1, tau2, tau3);

            // Check if SV is zero, use SH solution instead
            if (uSV.norm() < 1e-14) {
                ux = -py; uy = px; uz = 0;
                tau1 = c44 * py * qi;
                tau2 = -c44 * px * qi;
                tau3 = 0;
                uSV = Vector3c(ux, uy, uz);
                tauSV = Vector3c(tau1, tau2, tau3);
            }

            Complex flux = (tauSV.adjoint() * uSV).value() + (uSV.adjoint() * tauSV).value();
            Real epsilon_SV = 1.0 / std::sqrt(std::abs(flux) + 1e-30);

            eigvec.head<3>() = epsilon_SV * uSV;
            eigvec.tail<3>() = epsilon_SV * tauSV;

        } else {
            // SH-wave modes (i == 2 || i == 5)
            Real ux = -py, uy = px, uz = 0;
            Real tau1 = c44 * py * qi;
            Real tau2 = -c44 * px * qi;
            Real tau3 = 0;

            Vector3c uSH(ux, uy, uz);
            Vector3c tauSH(tau1, tau2, tau3);

            // Check if SH is zero, use SV solution instead
            if (uSH.norm() < 1e-14) {
                Real denom = (c13 + c44) * qi;
                if (std::abs(denom) < 1e-20) denom = 1e-20;
                Real L_SV = -(c11 * px*px + c44 * qi*qi - rho) / denom;

                ux = px; uy = py; uz = L_SV;
                tau1 = -c44 * px * (L_SV + qi);
                tau2 = -c44 * py * (L_SV + qi);
                tau3 = -(c13 * px*px + c33 * qi * L_SV);
                uSH = Vector3c(ux, uy, uz);
                tauSH = Vector3c(tau1, tau2, tau3);
            }

            Complex flux = (tauSH.adjoint() * uSH).value() + (uSH.adjoint() * tauSH).value();
            Real epsilon_SH = 1.0 / std::sqrt(std::abs(flux) + 1e-30);

            eigvec.head<3>() = epsilon_SH * uSH;
            eigvec.tail<3>() = epsilon_SH * tauSH;
        }

        result.D.col(i) = eigvec;
    }

    return result;
}

// ============================================================================
// compute_Rpp_aniso implementation (main algorithm)
// ============================================================================

ReflectionCoefficients compute_Rpp_aniso(
    Real freq,
    const std::vector<Real>& angles,
    const LayerModel& model,
    Real phi
) {
    const int nL = model.nLayers;
    const int Na = static_cast<int>(angles.size());

    ReflectionCoefficients result;
    result.Rpp.resize(Na);
    result.Rpsv.resize(Na);
    result.Rpsh.resize(Na);

    // Convert azimuth to radians
    Real phi_rad = deg2rad(phi);

    // Reference slowness (based on top layer vertical P velocity)
    Real vpz_top = std::sqrt(model.stiffness[0](2,2) / model.density[0]);
    Real Sp = 1.0 / vpz_top;

    // Angular frequency (real for phase calculation)
    Real w_real = 2.0 * PI * freq;

    // Precompute T0 (eigenvalues at normal incidence for phase calculation)
    std::vector<Matrix6c> T0(nL);
    for (int ii = 0; ii < nL; ++ii) {
        EigenResult eigres = compute_eig_VTI(0, 0, model.density[ii], model.stiffness[ii]);
        T0[ii] = eigres.T;
    }

    // Precompute safe slowness limit to avoid supercritical angles
    // MATLAB: pmax_safe = 0.95 * min(1./vpz) = 0.95 / max(vpz)
    Real vpz_max = 0.0;
    for (int ii = 0; ii < nL; ++ii) {
        Real vpz = std::sqrt(model.stiffness[ii](2,2) / model.density[ii]);
        vpz_max = std::max(vpz_max, vpz);
    }
    Real pmax_safe = 0.95 / vpz_max;

    // Main loop over incidence angles (OpenMP parallel)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Na; ++i) {
        Real theta = deg2rad(angles[i]);
        Real px = Sp * std::sin(theta) * std::cos(phi_rad);
        Real py = Sp * std::sin(theta) * std::sin(phi_rad);

        // Soft cap slowness to avoid supercritical angles
        Real p = std::hypot(px, py);
        if (p > pmax_safe && p > 0) {
            Real scale_p = pmax_safe / p;
            px *= scale_p;
            py *= scale_p;
        }

        // Storage for this angle
        std::vector<Matrix6c> D(nL), T(nL);
        std::vector<Matrix3c> Ed(nL-1), Eu(nL-1);
        std::vector<Matrix3c> ru(nL), rd(nL), tu(nL), td(nL);
        std::vector<Matrix3c> Ru(nL), Rd(nL), Tu(nL), Td(nL);

        // Initialize
        for (int ii = 0; ii < nL; ++ii) {
            ru[ii] = Matrix3c::Zero();
            rd[ii] = Matrix3c::Zero();
            tu[ii] = Matrix3c::Zero();
            td[ii] = Matrix3c::Zero();
            Ru[ii] = Matrix3c::Zero();
            Rd[ii] = Matrix3c::Zero();
            Tu[ii] = Matrix3c::Zero();
            Td[ii] = Matrix3c::Zero();
        }
        Tu[0] = Matrix3c::Identity();
        Td[0] = Matrix3c::Identity();

        // Compute eigenvalues/eigenvectors for each layer
        for (int ii = 0; ii < nL; ++ii) {
            EigenResult eigres = compute_eig_aniso(px, py, model.density[ii], model.stiffness[ii]);
            D[ii] = eigres.D;
            T[ii] = eigres.T;
        }

        // Compute phase shift matrices (Eu, Ed)
        for (int ii = 0; ii < nL-1; ++ii) {
            Real h = model.thickness[ii];
            Vector6c q0;
            for (int k = 0; k < 6; ++k) {
                q0(k) = T0[ii](k, k);
            }

            Vector6c e;
            for (int k = 0; k < 6; ++k) {
                e(k) = std::exp(q0(k) * h * w_real * Complex(0, 1));
            }

            Eu[ii] = Matrix3c::Zero();
            Ed[ii] = Matrix3c::Zero();
            for (int k = 0; k < 3; ++k) {
                Eu[ii](k, k) = e(k);
                Ed[ii](k, k) = e(k + 3);
            }
        }

        // Compute interface reflection/transmission matrices
        for (int ii = 1; ii < nL; ++ii) {
            // Q = D_below^(-1) * D_above
            Matrix6c Q = D[ii].fullPivLu().solve(D[ii-1]);

            Matrix3c Q11 = Q.block<3,3>(0, 0);
            Matrix3c Q12 = Q.block<3,3>(0, 3);
            Matrix3c Q21 = Q.block<3,3>(3, 0);
            Matrix3c Q22 = Q.block<3,3>(3, 3);

            Matrix3c invA11 = stable_inv(Q11);

            ru[ii] = Q21 * invA11;
            rd[ii] = -stable_solve(Q11, Q12);
            tu[ii] = invA11;
            td[ii] = Q22 - Q21 * invA11 * Q12;
        }

        // Recursive computation from surface to bottom
        Matrix3c I3 = Matrix3c::Identity();

        for (int ii = 1; ii < nL; ++ii) {
            int jj = ii - 1;

            // M1 = I - Eu^(-1) * rd * Ed * Ru
            Matrix3c Eu_inv = stable_inv(Eu[jj]);
            Matrix3c M1 = I3 - Eu_inv * rd[ii] * Ed[jj] * Ru[jj];
            Matrix3c M1i = stable_inv(M1);
            Tu[ii] = Tu[jj] * M1i * Eu_inv * tu[ii];

            // M2 = I - Ru * Eu^(-1) * rd * Ed
            Matrix3c M2 = I3 - Ru[jj] * Eu_inv * rd[ii] * Ed[jj];
            Matrix3c M2i = stable_inv(M2);
            Rd[ii] = Rd[jj] + Tu[jj] * Eu_inv * rd[ii] * Ed[jj] * M2i * Td[jj];

            // M3 = M1 (same)
            Matrix3c M3i = M1i;
            Ru[ii] = ru[ii] + td[ii] * Ed[jj] * Ru[jj] * M3i * Eu_inv * tu[ii];

            // M4 = M2 (same)
            Matrix3c M4i = M2i;
            Td[ii] = td[ii] * Ed[jj] * M4i * Td[jj];
        }

        // Extract reflection coefficients (conjugate to match MATLAB convention)
        result.Rpp[i] = std::conj(Rd[nL-1](0, 0));
        result.Rpsv[i] = std::conj(Rd[nL-1](0, 1));
        result.Rpsh[i] = std::conj(Rd[nL-1](0, 2));
    }

    return result;
}

// ============================================================================
// Batch computation implementation
// ============================================================================

std::vector<ReflectionCoefficients> compute_Rpp_aniso_batch(
    const std::vector<Real>& frequencies,
    const std::vector<Real>& angles,
    const LayerModel& model,
    Real phi
) {
    const int nFreq = static_cast<int>(frequencies.size());

    std::vector<ReflectionCoefficients> results(nFreq);

    // Parallel over frequencies
    #pragma omp parallel for schedule(dynamic)
    for (int f = 0; f < nFreq; ++f) {
        results[f] = compute_Rpp_aniso(frequencies[f], angles, model, phi);
    }

    return results;
}

// ============================================================================
// C interface implementation
// ============================================================================

extern "C" {

void azrm_compute_Rpp(
    double freq,
    const double* angles,
    int nAngles,
    const double* thickness,
    const double* density,
    const double* stiffness,
    int nLayers,
    double phi,
    double* Rpp_real,
    double* Rpp_imag,
    double* Rpsv_real,
    double* Rpsv_imag,
    double* Rpsh_real,
    double* Rpsh_imag
) {
    // Build LayerModel
    LayerModel model(nLayers);
    for (int i = 0; i < nLayers; ++i) {
        model.thickness[i] = thickness[i];
        model.density[i] = density[i];

        // Extract 6x6 stiffness matrix (column-major from MATLAB)
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                model.stiffness[i](row, col) = stiffness[row + col*6 + i*36];
            }
        }
    }

    // Build angles vector
    std::vector<Real> anglesVec(angles, angles + nAngles);

    // Compute
    ReflectionCoefficients result = compute_Rpp_aniso(freq, anglesVec, model, phi);

    // Copy results
    for (int i = 0; i < nAngles; ++i) {
        Rpp_real[i] = result.Rpp[i].real();
        Rpp_imag[i] = result.Rpp[i].imag();
        Rpsv_real[i] = result.Rpsv[i].real();
        Rpsv_imag[i] = result.Rpsv[i].imag();
        Rpsh_real[i] = result.Rpsh[i].real();
        Rpsh_imag[i] = result.Rpsh[i].imag();
    }
}

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
    double* Rpp_real,
    double* Rpp_imag,
    double* Rpsv_real,
    double* Rpsv_imag,
    double* Rpsh_real,
    double* Rpsh_imag
) {
    // Build LayerModel
    LayerModel model(nLayers);
    for (int i = 0; i < nLayers; ++i) {
        model.thickness[i] = thickness[i];
        model.density[i] = density[i];

        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                model.stiffness[i](row, col) = stiffness[row + col*6 + i*36];
            }
        }
    }

    std::vector<Real> anglesVec(angles, angles + nAngles);
    std::vector<Real> freqVec(frequencies, frequencies + nFreq);

    std::vector<ReflectionCoefficients> results = compute_Rpp_aniso_batch(freqVec, anglesVec, model, phi);

    // Copy results (row-major: [freq][angle])
    for (int f = 0; f < nFreq; ++f) {
        for (int i = 0; i < nAngles; ++i) {
            int idx = f * nAngles + i;
            Rpp_real[idx] = results[f].Rpp[i].real();
            Rpp_imag[idx] = results[f].Rpp[i].imag();
            Rpsv_real[idx] = results[f].Rpsv[i].real();
            Rpsv_imag[idx] = results[f].Rpsv[i].imag();
            Rpsh_real[idx] = results[f].Rpsh[i].real();
            Rpsh_imag[idx] = results[f].Rpsh[i].imag();
        }
    }
}

int azrm_get_num_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

void azrm_set_num_threads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
}

} // extern "C"

} // namespace azrm
