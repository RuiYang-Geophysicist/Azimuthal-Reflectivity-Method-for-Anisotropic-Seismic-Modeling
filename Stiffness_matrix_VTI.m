function [C_VTI] = Stiffness_matrix_VTI(model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the stiffness tensor C_ij for VTI media based on Thomsen parameters
%
% INPUT:
%   model : [N x 7] matrix with each row:
%           [thickness, vp, vs, rho, epsilon, delta, gamma]
%
% OUTPUT:
%   C_VTI : 6x6xN stiffness tensor for N layers (in Voigt notation)
%
% Reference:
%   Thomsen, L. (1986). Weak elastic anisotropy. Geophysics, 51(10), 1954–1966.
%
%   Written by Rui Yang, Tongji University, 2024/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract elastic and anisotropic parameters
vp      = model(:, 2);     % P-wave velocity (km/s)
vs      = model(:, 3);     % S-wave velocity (km/s)
rho     = model(:, 4);     % Density (g/cm^3)
epsilon = model(:, 5);     % Thomsen parameter
delta   = model(:, 6);     % Thomsen parameter
gamma   = model(:, 7);     % Thomsen parameter
nL      = size(model, 1);  % Number of layers

% Compute stiffness coefficients (isotropic + Thomsen)
C33 = vp.^2 .* rho;
C44 = vs.^2 .* rho;
C55 = C44;

C11 = C33 .* (1 + 2 .* epsilon);
C13 = sqrt(2 .* delta .* C33 .* (C33 - C55) + (C33 - C55).^2) - C55;
C12 = C11 - 2 .* C55;
C66 = C44 .* (1 + 2 .* gamma);

% Initialize stiffness tensor
C_VTI = zeros(6, 6, nL);

% Construct 6x6 symmetric stiffness matrix for each layer
for c = 1:nL
    C_VTI(:, :, c) = [ C11(c), C12(c), C13(c),     0,     0,     0;
                       C12(c), C11(c), C13(c),     0,     0,     0;
                       C13(c), C13(c), C33(c),     0,     0,     0;
                            0,      0,      0, C44(c),     0,     0;
                            0,      0,      0,     0, C44(c),     0;
                            0,      0,      0,     0,     0, C66(c)];
end

end
