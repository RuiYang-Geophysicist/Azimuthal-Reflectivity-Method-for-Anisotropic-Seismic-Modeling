function [Rpp, Rpsv, Rpsh] = compute_Rpp_aniso(f, angle, thickness, c_layer, den, phi)
% COMPUTE_RPP_ANISO Compute PP, PSv, PSh reflection coefficients for general anisotropic media
%
% This function uses the reflectivity method (Fryer & Frazer, 1984, 1987) with
% numerical eigenvalue decomposition to compute reflection coefficients for
% arbitrary anisotropic layered media.
%
% INPUTS:
%   f         : Frequency (Hz)
%   angle     : Incidence angles (degrees), array
%   thickness : Layer thicknesses (km), array [nL x 1]
%   c_layer   : Stiffness tensors (6x6xnL)
%   den       : Densities (g/cm^3), array [nL x 1]
%   phi       : Azimuth angle (degrees)
%
% OUTPUTS:
%   Rpp  : PP reflection coefficients at each angle
%   Rpsv : P-to-SV reflection coefficients
%   Rpsh : P-to-SH reflection coefficients
%
% SUPPORTED MEDIA:
%   - Isotropic
%   - VTI (Vertical Transverse Isotropy)
%   - HTI (Horizontal Transverse Isotropy)
%   - Orthorhombic
%   - Monoclinic
%   - Triclinic (fully general)
%
% REFERENCES:
%   Fryer & Frazer, GJI 1984, 1987
%   Kennett, 1983
%
% Written by Rui Yang, Tongji University, 2024
% General anisotropy version using numerical eigenvalue solver

%% Initialize
nL = length(den);              % Number of layers
Na = length(angle);            % Number of incidence angles
phi = phi * pi / 180;          % Convert azimuth to radians

% Reference slowness (based on top layer vertical P velocity)
Sp = 1 / sqrt(c_layer(3,3,1) / den(1));

%% Preallocate matrices
% Eigenvalue/eigenvector matrices
T0 = zeros(6, 6, nL);          % Eigenvalues at normal incidence (for phase)
D  = zeros(6, 6, nL);          % Eigenvectors at current angle
T  = zeros(6, 6, nL);          % Eigenvalues at current angle

% Phase shift matrices
Ed = zeros(3, 3, nL-1);
Eu = zeros(3, 3, nL-1);

% Interface reflection/transmission matrices
ru = zeros(3, 3, nL);
rd = zeros(3, 3, nL);
tu = zeros(3, 3, nL);
td = zeros(3, 3, nL);

% Cumulative reflection/transmission matrices
Ru = zeros(3, 3, nL);
Rd = zeros(3, 3, nL);
Tu = zeros(3, 3, nL);
Td = zeros(3, 3, nL);
Tu(:,:,1) = eye(3);
Td(:,:,1) = eye(3);

% Q matrix blocks
Q   = zeros(6, 6, nL);
Q11 = zeros(3, 3, nL);
Q12 = zeros(3, 3, nL);
Q21 = zeros(3, 3, nL);
Q22 = zeros(3, 3, nL);

% Output arrays
Rpp  = zeros(1, Na);
Rpsv = zeros(1, Na);
Rpsh = zeros(1, Na);

%% Angular frequency
% Use real frequency for phase calculation to preserve amplitude
w_real = 2 * pi * f;

%% Compute T0 (eigenvalues at normal incidence for phase calculation)
% This gives "flat gathers" similar to NMO-corrected data
% IMPORTANT: Use VTI analytical solver for T0 to ensure consistent ordering
% At normal incidence, SV and SH are degenerate (same eigenvalue), and
% numerical solver may order them differently. VTI solver handles this correctly.
px = 0;
py = 0;
for ii = 1:nL
    % Use VTI solver for normal incidence (handles degeneracy correctly)
    [~, T0(:,:,ii)] = compute_eig_VTI(px, py, den(ii), c_layer(:,:,ii));
end

%% Precompute safe slowness limit to avoid supercritical angles
vpz = zeros(1, nL);
for ii = 1:nL
    vpz(ii) = sqrt(c_layer(3,3,ii) / den(ii));
end
pmax_safe = 0.95 * min(1 ./ vpz);

%% Main loop over incidence angles
for i = 1:Na
    theta = angle(i) * pi / 180;
    px = Sp * sin(theta) * cos(phi);
    py = Sp * sin(theta) * sin(phi);

    % Soft cap slowness to avoid supercritical angles
    p = hypot(px, py);
    if p > pmax_safe && p > 0
        scale_p = pmax_safe / p;
        px = px * scale_p;
        py = py * scale_p;
    end

    %% Compute eigenvalues/eigenvectors for each layer (numerical method)
    for ii = 1:nL
        [D(:,:,ii), T(:,:,ii)] = compute_eig_aniso(px, py, den(ii), c_layer(:,:,ii));
    end

    %% Compute phase shift matrices (Eu, Ed)
    % Use T0 (normal incidence) for flat gathers
    % Use real frequency to preserve amplitude
    for ii = 1:nL-1
        h = thickness(ii);
        q0 = diag(T0(:,:,ii));
        e = exp(q0 * h * w_real * 1i);
        Eu(:,:,ii) = diag(e(1:3));
        Ed(:,:,ii) = diag(e(4:6));
    end

    %% Compute interface reflection/transmission matrices
    for ii = 2:nL
        % Q = D_below^(-1) * D_above (Fryer & Frazer 1984, Eq. 4.22)
        Q(:,:,ii) = D(:,:,ii) \ D(:,:,ii-1);

        Q11(:,:,ii) = Q(1:3, 1:3, ii);
        Q12(:,:,ii) = Q(1:3, 4:6, ii);
        Q21(:,:,ii) = Q(4:6, 1:3, ii);
        Q22(:,:,ii) = Q(4:6, 4:6, ii);

        % Stable matrix operations
        A11 = Q11(:,:,ii);
        invA11 = stable_inv(A11);

        ru(:,:,ii) = Q21(:,:,ii) * invA11;
        rd(:,:,ii) = -stable_solve(A11, Q12(:,:,ii));
        tu(:,:,ii) = invA11;
        td(:,:,ii) = Q22(:,:,ii) - Q21(:,:,ii) * invA11 * Q12(:,:,ii);
    end

    %% Recursive computation from surface to bottom
    I3 = eye(3);
    for ii = 2:nL
        jj = ii - 1;

        % M1 = I - Eu^(-1) * rd * Ed * Ru
        M1 = I3 - (Eu(:,:,jj) \ rd(:,:,ii)) * Ed(:,:,jj) * Ru(:,:,jj);
        M1i = stable_inv(M1);
        Tu(:,:,ii) = Tu(:,:,jj) * M1i / Eu(:,:,jj) * tu(:,:,ii);

        % M2 = I - Ru * Eu^(-1) * rd * Ed
        M2 = I3 - (Ru(:,:,jj) / Eu(:,:,jj)) * rd(:,:,ii) * Ed(:,:,jj);
        M2i = stable_inv(M2);
        Rd(:,:,ii) = Rd(:,:,jj) + Tu(:,:,jj) / Eu(:,:,jj) * rd(:,:,ii) * Ed(:,:,jj) * M2i * Td(:,:,jj);

        % M3 = M1 (same)
        M3i = M1i;
        Ru(:,:,ii) = ru(:,:,ii) + td(:,:,ii) * Ed(:,:,jj) * Ru(:,:,jj) * M3i / Eu(:,:,jj) * tu(:,:,ii);

        % M4 = M2 (same)
        M4i = M2i;
        Td(:,:,ii) = td(:,:,ii) * Ed(:,:,jj) * M4i * Td(:,:,jj);
    end

    %% Extract reflection coefficients
    % Note: Use conjugate transpose (') to match compute_Rpp_VTI convention
    Rpp(i)  = Rd(1, 1, nL)';
    Rpsv(i) = Rd(1, 2, nL)';
    Rpsh(i) = Rd(1, 3, nL)';
end

end

%% ========== Local Helper Functions ==========

function X = stable_solve(A, B)
% STABLE_SOLVE Solve A*X = B with regularization for near-singular A

    % Regularize if needed
    A = regularize_matrix(A);

    % Check for invalid values
    if any(~isfinite(A(:)))
        A(~isfinite(A)) = 0;
        A = A + 1e-6 * eye(size(A));
    end

    % Solve
    X = A \ B;

    % Fallback to pinv if result is invalid
    if any(~isfinite(X(:)))
        X = pinv(A) * B;
    end
end

function X = stable_inv(A)
% STABLE_INV Compute inverse with regularization for near-singular A

    % Regularize if needed
    A = regularize_matrix(A);

    % Check for invalid values
    if any(~isfinite(A(:)))
        A(~isfinite(A)) = 0;
        A = A + 1e-6 * eye(size(A));
    end

    % Try LU decomposition first
    [L, U, P] = lu(A);

    if any(~isfinite(L(:))) || any(~isfinite(U(:))) || any(abs(diag(U)) < 1e-14)
        % Fallback to pseudoinverse
        X = pinv(A);
    else
        X = U \ (L \ P);
    end
end

function A = regularize_matrix(A)
% REGULARIZE_MATRIX Add small regularization if matrix is near-singular

    rc = rcond(A);
    if ~isfinite(rc) || rc < 1e-12
        A = A + 1e-6 * norm(A, 'fro') * eye(size(A));
    end
end
