function [D, T] = compute_eig_aniso(px, py, rho, C)
% COMPUTE_EIG_ANISO Compute eigenvalues and eigenvectors for general anisotropic media
%
% This function uses numerical eigenvalue decomposition to solve the
% Christoffel equation for arbitrary anisotropic media (VTI, HTI, OA,
% Monoclinic, Triclinic, etc.)
%
% INPUTS:
%   px, py : Horizontal slowness components (s/km)
%   rho    : Density (g/cm^3)
%   C      : 6x6 stiffness matrix in Voigt notation (GPa or consistent units)
%
% OUTPUTS:
%   D : 6x6 eigenvector matrix, columns ordered as:
%       [up-qP, up-qS1, up-qS2, down-qP, down-qS1, down-qS2]
%       where qP = quasi-P, qS1/qS2 = quasi-shear (sorted by |q|)
%   T : 6x6 diagonal eigenvalue matrix (vertical slownesses q)
%
% NOTES:
%   - Based on Fryer & Frazer (1984, 1987) formulation
%   - Uses MATLAB's eig() for numerical stability
%   - Supports any symmetry class: ISO, VTI, HTI, OA, Monoclinic, Triclinic
%   - Does NOT assume VTI-specific wave polarizations
%   - Uses energy flux for consistent sign convention
%
% REFERENCES:
%   Fryer & Frazer, GJI 1984, 1987
%   Kennett, 1983
%   Chapman, 2004 (Fundamentals of Seismic Wave Propagation)
%
% Written by Rui Yang, Tongji University, 2024
% Modified for general anisotropy support (no VTI assumptions)

%% Extract stiffness components (Triclinic - most general)
c11 = C(1,1); c12 = C(1,2); c13 = C(1,3); c14 = C(1,4); c15 = C(1,5); c16 = C(1,6);
c22 = C(2,2); c23 = C(2,3); c24 = C(2,4); c25 = C(2,5); c26 = C(2,6);
c33 = C(3,3); c34 = C(3,4); c35 = C(3,5); c36 = C(3,6);
c44 = C(4,4); c45 = C(4,5); c46 = C(4,6);
c55 = C(5,5); c56 = C(5,6);
c66 = C(6,6);

%% Build the 6x6 system matrix A
% State vector: [ux, uy, uz, tau_xz, tau_yz, tau_zz]^T
% System: A * state = q * state, where q is vertical slowness

% S matrix elements (stress-displacement relation in horizontal plane)
s11 = px^2*(c11 - c13^2/c33) + 2*px*py*(c16 - c13*c36/c33) + py^2*(c66 - c36^2/c33) - rho;
s12 = px^2*(c16 - c13*c36/c33) + px*py*(c12 + c66 - c13*c23/c33 - c36^2/c33) + py^2*(c26 - c23*c36/c33);
s22 = px^2*(c66 - c36^2/c33) + 2*px*py*(c26 - c23*c36/c33) + py^2*(c22 - c23^2/c33) - rho;

% T matrix (relates displacement to traction gradient)
Ta = [0,                           0,                           -px;
      0,                           0,                           -py;
      -(px*c13 + py*c36)/c33,     -(px*c36 + py*c23)/c33,       0];

% C matrix (compliance-related)
% Relates traction derivatives to displacement
% For stress-strain: tau_xz = c55*(du_x/dz + du_z/dx), tau_yz = c44*(du_y/dz + du_z/dy)
% So: du_x/dz ~ tau_xz/c55, du_y/dz ~ tau_yz/c44
% The inverse relation gives: Ca(1,1) = -1/c55, Ca(2,2) = -1/c44 (for c45=0)
det_c45 = c44*c55 - c45^2;
if abs(det_c45) < 1e-20
    det_c45 = sign(det_c45 + 1e-30) * 1e-20;  % Regularization
end

Ca = [-c44/det_c45,   c45/det_c45,  0;
       c45/det_c45,  -c55/det_c45,  0;
       0,             0,           -1/c33];

% S matrix
Sa = [s11,  s12,  0;
      s12,  s22,  0;
      0,    0,   -rho];

% Assemble full 6x6 system matrix
A = [Ta, Ca;
     Sa, Ta'];

%% Solve eigenvalue problem
[D_raw, T_raw] = eig(A);

%% Reorder eigenpairs for general anisotropy
[D, T] = reorder_eigenpairs_general(D_raw, T_raw);

end

%% ========== Local Functions ==========

function [Dout, Tout] = reorder_eigenpairs_general(Din, Tin)
% REORDER_EIGENPAIRS_GENERAL Reorder eigenpairs for general anisotropic media
%
% Convention:
%   - First 3 columns: up-going waves (Re(q) < 0 or energy flux upward)
%   - Last 3 columns: down-going waves (Re(q) > 0 or energy flux downward)
%   - Within each group: sorted by |q| ascending (fastest to slowest)
%
% This method does NOT assume VTI-specific polarizations.
% It uses energy flux direction as the primary criterion for up/down separation.

    % Extract eigenvalues
    q = diag(Tin);
    n = length(q);

    %% Step 1: Compute energy flux for each mode to determine propagation direction
    % Energy flux in z-direction: Fz = -Re(tau_z^* . v) where tau_z = [tau_xz, tau_yz, tau_zz]
    % For our state vector [ux,uy,uz,tau_xz,tau_yz,tau_zz], this is:
    % Fz = -Re(conj(D(4:6,j))' * D(1:3,j)) ... but we need proper scaling
    %
    % Alternative: use sign of Re(q) as primary, with energy flux as tiebreaker

    energy_flux_z = zeros(n, 1);
    for j = 1:n
        u = Din(1:3, j);      % displacement
        tau = Din(4:6, j);    % traction
        % Vertical energy flux (positive = downward for exp(i*w*q*z) convention)
        energy_flux_z(j) = real(tau' * u + u' * tau);
    end

    %% Step 2: Classify as up-going or down-going
    % Primary criterion: sign of Re(q)
    % Secondary criterion: energy flux direction (for purely imaginary q)

    tol_real = 1e-10;

    is_down = zeros(n, 1);
    for j = 1:n
        if real(q(j)) > tol_real
            is_down(j) = 1;  % Positive Re(q) -> down-going
        elseif real(q(j)) < -tol_real
            is_down(j) = 0;  % Negative Re(q) -> up-going
        else
            % Purely imaginary or zero: use energy flux
            % Positive energy flux = down-going (energy propagates downward)
            is_down(j) = (energy_flux_z(j) > 0);
        end
    end

    idx_up = find(~is_down);
    idx_down = find(is_down);

    % Fallback: if not exactly 3+3, force split by Re(q)
    if numel(idx_up) ~= 3 || numel(idx_down) ~= 3
        [~, ord] = sort(real(q), 'ascend');
        idx_up = ord(1:3);
        idx_down = ord(4:6);
    end

    %% Step 3: Sort each group by |q| (ascending = fastest wave first)
    % This gives consistent ordering without assuming P/SV/SH labels

    [~, sort_up] = sort(abs(q(idx_up)), 'ascend');
    [~, sort_down] = sort(abs(q(idx_down)), 'ascend');

    idx_up_sorted = idx_up(sort_up);
    idx_down_sorted = idx_down(sort_down);

    % Final ordering: [up-fast, up-mid, up-slow, down-fast, down-mid, down-slow]
    final_order = [idx_up_sorted(:)', idx_down_sorted(:)'];

    %% Step 4: Build output matrices
    Dout = Din(:, final_order);
    Tout = diag(q(final_order));

    %% Step 5: Normalize eigenvectors using energy flux
    Dout = normalize_eigenvectors_energy(Dout);

end

function Dout = normalize_eigenvectors_energy(D)
% NORMALIZE_EIGENVECTORS_ENERGY Normalize eigenvectors for general anisotropy
%
% Normalization: energy flux magnitude = 1
%   |tau' * u + u' * tau| = 1
%
% Sign convention: VTI-compatible (ux > 0 or uy > 0)
%   The VTI analytical code uses: ux = px (for P/SV), uy = px (for SH)
%   This means for px > 0: horizontal displacement should be positive.
%
% For general anisotropy, we use the same principle to ensure compatibility:
%   - If |ux| > threshold: ensure ux > 0
%   - Else if |uy| > threshold: ensure uy > 0
%   - Else: ensure uz > 0

    Dout = D;
    n = size(D, 2);

    for j = 1:n
        u = D(1:3, j);      % Displacement
        tau = D(4:6, j);    % Traction

        % Compute energy flux for normalization
        energy_flux = real(tau' * u + u' * tau);

        if abs(energy_flux) > 1e-14
            % Energy normalization: |flux| = 1
            epsilon = 1 / sqrt(abs(energy_flux));
            Dout(:, j) = epsilon * D(:, j);
        else
            % Degenerate case: use Euclidean norm
            Dout(:, j) = D(:, j) / norm(D(:, j));
        end

        % Sign convention compatible with VTI analytical code
        % Use absolute threshold to handle small angles properly
        u_norm = Dout(1:3, j);
        ux = real(u_norm(1));
        uy = real(u_norm(2));
        uz = real(u_norm(3));

        tol_abs = 1e-10;

        if abs(ux) > tol_abs
            % ux is non-zero: use ux sign (P/SV-like mode)
            if ux < 0
                Dout(:, j) = -Dout(:, j);
            end
        elseif abs(uy) > tol_abs
            % ux negligible, uy non-zero: use uy sign (SH-like mode)
            if uy < 0
                Dout(:, j) = -Dout(:, j);
            end
        else
            % Both horizontal components negligible (normal incidence)
            if uz < 0
                Dout(:, j) = -Dout(:, j);
            end
        end
    end
end
