function [Npsi_1, Nphi_1] = Majorana_rebuild(psi_1,N_tot)
% Orthonormalize the degenerate zero-energy states.
% psi_1 is the wave-function for the first energy level.
%% basis
ty = [0 -1i; 1i 0];
syty = kron(ty,ty);
main_diag_vec = [ones(1,N_tot)];
main_diag = diag(main_diag_vec);
PHS = kron(main_diag,syty); % particle-hole symmetry operator

%% Rebuild first-level states
phi_1 = PHS*conj(psi_1); % particle-hole conjugate wave-function
gamma_1 = psi_1+phi_1;
norm1 = sqrt(gamma_1'*gamma_1);
gamma_1 = gamma_1./norm1;
gamma_2 = -1i.*(psi_1-phi_1);
norm2 = sqrt(gamma_2'*gamma_2);
gamma_2 = gamma_2./norm2;
Npsi_1 = (gamma_1+1i.*gamma_2)./sqrt(2); % Quasi-particle wave-function
Nphi_1 = (gamma_1-1i.*gamma_2)./sqrt(2); % Quasi-hole wave-function
