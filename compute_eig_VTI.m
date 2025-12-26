function [D, T] = compute_eig_VTI(px, py, rho, c)
% compute_eig_VTI computes the eigenvalue and eigenvector matrices of each
% layer using the analytical solutions by Fryer and Frazer (1987)
%
%   Written by Rui Yang, Tongji University, 2024/06
%              Neng Lu, Jilin University, 2015/05

% Initialize output
D = zeros(6, 6);
T = zeros(6, 6);

% Extract stiffness coefficients for VTI medium
c11 = c(1,1); c12 = c(1,2); c13 = c(1,3);
c33 = c(3,3); c44 = c(4,4); c66 = c(6,6);

% Auxiliary scalar terms (equation terms)
s11 = px^2 * (c11 - c13^2 / c33) + py^2 * c66 - rho;
s12 = px * py * (c12 + c66 - c13^2 / c33);
s22 = px^2 * c66 + py^2 * (c11 - c13^2 / c33) - rho;

% Construct system matrix A (Fryer & Frazer 1987, Eq. 4.9)
A1 = [0 0 -px*c13/c33 s11 s12 0]';  
A2 = [0 0 -py*c13/c33 s12 s22 0]';
A3 = [-px -py 0 0 0 -rho]';
A4 = [-1/c44 0 0 0 0 -px]';
A5 = [0 -1/c44 0 0 0 -py]';
A6 = [0 0 -1/c33 -px*c13/c33 -py*c13/c33 0]';
A  = [A1 A2 A3 A4 A5 A6];

% Solve eigenproblem
[~, eigvalues] = eig(A);
qi = real(diag(eigvalues));

% Sort eigenvalues: [qD_T > qD_S > qD_P > 0 > qU_P > qU_S > qU_T]
sorted_qi = sort(qi, 'descend');
q = [sorted_qi(4); sorted_qi(5); sorted_qi(6); sorted_qi(3); sorted_qi(2); sorted_qi(1)];

% Compute eigenvectors for each eigenvalue
b = zeros(6, 6);  % Preallocate

for i = 1:6
    if i == 1 || i == 4  % P-wave modes
        L_P  = -(c11 * px^2 + c44 * q(i)^2 - rho) / ((c13 + c44) * q(i));
        ux   = px; uy = py; uz = L_P;
        tau1 = -c44 * px * (L_P + q(i));
        tau2 = -c44 * py * (L_P + q(i));
        tau3 = -(c13 * px^2 + c33 * q(i) * L_P);
        
        uP   = [ux; uy; uz];
        tauP = [tau1; tau2; tau3];
        epsilon_p = 1 / sqrt(abs(tauP' * uP + uP' * tauP));
        b(:, i) = epsilon_p * [uP; tauP];

    elseif i == 2 || i == 5  % SV-wave modes
        L_SV = -(c11 * px^2 + c44 * q(i)^2 - rho) / ((c13 + c44) * q(i));
        ux   = px; uy = py; uz = L_SV;
        tau1 = -c44 * px * (L_SV + q(i));
        tau2 = -c44 * py * (L_SV + q(i));
        tau3 = -(c13 * px^2 + c33 * q(i) * L_SV);
        
        uSV   = [ux; uy; uz];
        tauSV = [tau1; tau2; tau3];
        epsilon_SV = 1 / sqrt(abs(tauSV' * uSV + uSV' * tauSV));
        b(:, i) = epsilon_SV * [uSV; tauSV];

        if norm(uSV) == 0  % Use SH solution if SV is zero
            ux = -py; uy = px; uz = 0;
            tau1 =  c44 * py * q(i);
            tau2 = -c44 * px * q(i);
            tau3 = 0;
            uSV   = [ux; uy; uz];
            tauSV = [tau1; tau2; tau3];
            epsilon_SV = 1 / sqrt(abs(tauSV' * uSV + uSV' * tauSV));
            b(:, i) = epsilon_SV * [uSV; tauSV];
        end

    elseif i == 3 || i == 6  % SH-wave modes
        ux = -py; uy = px; uz = 0;
        tau1 =  c44 * py * q(i);
        tau2 = -c44 * px * q(i);
        tau3 = 0;
        
        uSH   = [ux; uy; uz];
        tauSH = [tau1; tau2; tau3];
        epsilon_SH = 1 / sqrt(abs(tauSH' * uSH + uSH' * tauSH));
        b(:, i) = epsilon_SH * [uSH; tauSH];

        if norm(uSH) == 0  % Use SV solution if SH is zero
            L_SV = -(c11 * px^2 + c44 * q(i)^2 - rho) / ((c13 + c44) * q(i));
            ux   = px; uy = py; uz = L_SV;
            tau1 = -c44 * px * (L_SV + q(i));
            tau2 = -c44 * py * (L_SV + q(i));
            tau3 = -(c13 * px^2 + c33 * q(i) * L_SV);
            uSH   = [ux; uy; uz];
            tauSH = [tau1; tau2; tau3];
            epsilon_SH = 1 / sqrt(abs(tauSH' * uSH + uSH' * tauSH));
            b(:, i) = epsilon_SH * [uSH; tauSH];
        end
    end
end

% Assign output
D = b;
T = diag(q);

end
