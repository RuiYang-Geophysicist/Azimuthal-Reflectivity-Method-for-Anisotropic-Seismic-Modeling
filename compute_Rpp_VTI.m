function [Rpp, Rpsv, Rpsh] = compute_Rpp_VTI(f, angle, thickness, c_layer, den, phi)
% compute_Rpp_VTI computes PP, PSv, PSh reflection coefficients using
% Fryer and Frazer's eigenstructure-based method in VTI media.
%
%   Written by Rui Yang, Tongji University, 2024/06
%              Neng Lu, Jilin University, 2015/05

nL   = length(den);              % Number of layers
Na   = length(angle);            % Number of incidence angles
phi  = phi * pi / 180;           % Convert to radians

% Optional method switch via global variable: 'default' or 'zref'
% If 'zref', use z_Reflectivity_P_incidence3 to compute Rpp directly
global RPP_METHOD
if ~isempty(which('z_Reflectivity_P_incidence3')) && exist('RPP_METHOD','var') && ischar(RPP_METHOD) && strcmpi(RPP_METHOD,'zref')
    try
        h = thickness(1:end-1);                 % use top (nL-1) thickness values
        omega_deg = phi * 180 / pi;             % convert back to degrees for z_Reflectivity
        [Rpp, Rpsv, Rpsh, ~] = z_Reflectivity_P_incidence3(f, angle, h, c_layer, den, omega_deg);
        return;
    catch
        % fall through to default method on any error
    end
end

% Horizontal slowness (based on top layer vertical P velocity)
Sp = 1 / sqrt(c_layer(3, 3, 1) / den(1));

% Preallocate matrices
T0(:,:,:) = zeros(6,6,nL);  % initialize eigenvalue T0 and eigenvector D0 matrices 
D0(:,:,:) = zeros(6,6,nL);  % as waves propagating vertically
T(:,:,:)  = zeros(6,6,nL);  % initialize eigenvalue T and eigenvector D matrices
D(:,:,:)  = zeros(6,6,nL);  % as waves propagating arbitrarily
Ed(:,:,:) = zeros(3,3,nL-1);% initialize the phase shift matrix Ed,Eu
Eu(:,:,:) = zeros(3,3,nL-1);
ru(:,:,:) = zeros(3,3,nL);  % initialize the reflection matrices ru,rd of i layer
rd(:,:,:) = zeros(3,3,nL);
tu(:,:,:) = zeros(3,3,nL);  % initialize the transmission matrices tu,td of i layer
td(:,:,:) = zeros(3,3,nL);
Ru(:,:,:) = zeros(3,3,nL);  % initialize the reflection matrices Ru,Rd of i+1 layer
Rd(:,:,:) = zeros(3,3,nL);  
Tu(:,:,1) = eye(3);                % initialize the transmission matrices Tu,Rd of i+1 layer
Td(:,:,1) = eye(3);    

Q    = zeros(6, 6, nL);
Q11  = zeros(3, 3, nL); Q12 = zeros(3, 3, nL);
Q21  = zeros(3, 3, nL); Q22 = zeros(3, 3, nL);

Rpp  = zeros(1, Na);
Rpsv = zeros(1, Na);
Rpsh = zeros(1, Na);

% Angular frequency with tiny complex damping (stabilizes resonance/evanescence)
% NOTE: Use separate real/complex frequencies for phase calculation
%       - w_real: for phase amplitude (|exp(i*w*q*h)| = 1)
%       - w_complex: for numerical stability in matrix operations
alpha = 0.5;            % Hz; damping for stability
w_real = 2 * pi * f;    % Real frequency for phase (preserves amplitude)
w = 2 * pi * (f + 1i*alpha);  % Complex frequency for stability

%% Compute T0 (eigenvalues at normal incidence)
px = 0;
py = 0;
for ii = 1:nL
    [~, T0(:, :, ii)] = compute_eig_VTI(px, py, den(ii), c_layer(:, :, ii));
end


% Precompute vertical P velocity per layer & global safe p-limit
vpz = zeros(1,nL);
for ii = 1:nL
    vpz(ii) = sqrt(c_layer(3,3,ii) / den(ii));   % km/s
end
pmax_safe = 0.95 * min(1 ./ vpz);  % s/km，全层最小的1/vp，乘0.95作软限幅

%% Loop over incidence angles
for i = 1:Na
    theta = angle(i) * pi / 180;
    px = Sp * sin(theta) * cos(phi);
    py = Sp * sin(theta) * sin(phi);
    
    % Soft cap p to avoid supercritical across all layers (keeps gathers horizontal)
    p = hypot(px, py);
    if p > pmax_safe && p > 0
        scale_p = pmax_safe / p;   % 0<scale_p<=1
        px = px * scale_p;
        py = py * scale_p;
    end

    % Compute eigenvalues/vectors for each layer
    for ii = 1:nL
        [D(:, :, ii), T(:, :, ii)] = compute_eig_VTI(px, py, den(ii), c_layer(:, :, ii));
    end

    %% Compute phase shift matrices (Eu, Ed) - vectorized inner loop
    % Use w_real for phase to preserve amplitude (|exp(i*w_real*q*h)| = 1)
    for ii = 1:nL-1
        h = thickness(ii);
        % Vectorized: extract diagonal and compute exp in one step
        e = exp(diag(T0(:, :, ii)) * h * w_real * 1i);
        Eu(:, :, ii) = diag(e(1:3));
        Ed(:, :, ii) = diag(e(4:6));
    end


    %% Compute reflection/transmission at each interface
    for ii = 2:nL
        Q(:, :, ii)   = D(:, :, ii) \ D(:, :, ii-1);  % Fryer & Frazer (1984), Eq. (4.22)
        Q11(:, :, ii) = Q(1:3, 1:3, ii);
        Q12(:, :, ii) = Q(1:3, 4:6, ii);
        Q21(:, :, ii) = Q(4:6, 1:3, ii);
        Q22(:, :, ii) = Q(4:6, 4:6, ii);

        % Stable solves for interface matrices (regularized)
        A11 = Q11(:, :, ii);
        invA11 = inv_stable3(A11);
        ru(:, :, ii) = Q21(:, :, ii) * invA11;                    % = Q21 / A11
        rd(:, :, ii) = - solve_stable3(A11, Q12(:, :, ii));       % = - A11 \\ Q12
        tu(:, :, ii) = invA11;
        td(:, :, ii) = Q22(:, :, ii) - Q21(:, :, ii) * invA11 * Q12(:, :, ii);

    end

    %% Recursive computation from free surface to last layer
    for ii = 2:nL
        jj = ii - 1;

    % Stable inverses for recursion blocks
    M1 = eye(3) - (Eu(:, :, jj) \ rd(:, :, ii)) * Ed(:, :, jj) * Ru(:, :, jj);
    M1i = inv_stable3(M1);
    Tu(:, :, ii) = Tu(:, :, jj) * M1i / Eu(:, :, jj) * tu(:, :, ii);
    
    M2 = eye(3) - (Ru(:, :, jj) / Eu(:, :, jj)) * rd(:, :, ii) * Ed(:, :, jj);
    M2i = inv_stable3(M2);
    Rd(:, :, ii) = Rd(:, :, jj) + Tu(:, :, jj) / Eu(:, :, jj) * rd(:, :, ii) * Ed(:, :, jj) * M2i * Td(:, :, jj);
    
    M3 = eye(3) - (Eu(:, :, jj) \ rd(:, :, ii)) * Ed(:, :, jj) * Ru(:, :, jj);
    M3i = inv_stable3(M3);
    Ru(:, :, ii) = ru(:, :, ii) + td(:, :, ii) * Ed(:, :, jj) * Ru(:, :, jj) * M3i / Eu(:, :, jj) * tu(:, :, ii);
    
    M4 = eye(3) - (Ru(:, :, jj) / Eu(:, :, jj)) * rd(:, :, ii) * Ed(:, :, jj);
    M4i = inv_stable3(M4);
    Td(:, :, ii) = td(:, :, ii) * Ed(:, :, jj) * M4i * Td(:, :, jj);

    end

    %% Extract reflection coefficients at surface
    Rpp(i)  = Rd(1, 1, nL)';
    Rpsv(i) = Rd(1, 2, nL)';
    Rpsh(i) = Rd(1, 3, nL)';
end

end

function X = solve_stable3(A, B)
    A = regularize3(A);
    if any(~isfinite(A(:))) || any(isnan(A(:)))
        A(~isfinite(A)) = 0;
        A = A + 1e-6*eye(3);
    end
    X = A \ B;
    if any(~isfinite(X(:))) || any(isnan(X(:)))
        X = pinv(A) * B;
    end
end

function X = inv_stable3(A)
    A = regularize3(A);
    if any(~isfinite(A(:))) || any(isnan(A(:)))
        A(~isfinite(A)) = 0;
        A = A + 1e-6*eye(3);
    end
    % Prefer backslash to explicit inv; fallback to pinv on failure
    [L,U,P] = lu(A);
    if any(~isfinite(L(:))) || any(~isfinite(U(:))) || any(isnan(U(:))) || any(abs(diag(U)) < 1e-14)
        X = pinv(A);
    else
        X = U \ (L \ P);
    end
end

function A = regularize3(A)
    rc = rcond(A);
    if ~isfinite(rc) || rc < 1e-10
        A = A + 3e-6 * norm(A,'fro') * eye(3);  % 可调：1e-6→3e-6 or 1e-5
    end
end