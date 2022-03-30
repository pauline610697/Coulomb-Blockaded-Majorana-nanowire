function H = hse_v10(t,Delta,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,omega)
%% Spin & Particle-hole spaces (Pauli matrices)
s0 = [1 0; 0 1]; sx = [0 1; 1 0]; sy = [0 -1i; 1i 0]; sz = [1 0; 0 -1]; % Spin space
t0 = [1 0; 0 1]; tx = [0 1; 1 0]; ty = [0 -1i; 1i 0]; tz = [1 0; 0 -1]; % particle-hole space

tzs0 = kron(tz,s0); tzsy = kron(tz,sy); t0s0 = kron(t0,s0); txs0 = kron(tx,s0); t0sx = kron(t0,sx);

band11sm = spdiags([ones(N_tot,1) ones(N_tot,1)],[-1,1],N_tot,N_tot);
band1m1sm = spdiags([ones(N_tot,1) -ones(N_tot,1)],[-1,1],N_tot,N_tot);
eyesm = speye(N_tot);

barrier = spdiags(cat(1,ones(Nbarrier,1),zeros(N_tot-Nbarrier,1)),0,N_tot,N_tot);

dot_vec1 = (Nbarrier+1):N_dot;
dotV1 = VD1*cos(1.5.*pi.*dot_vec1./N_dot)';
dot1 = spdiags(cat(1,zeros(Nbarrier,1),dotV1,zeros(N_tot-N_dot,1)),0,N_tot,N_tot);

dot_vec2 = 1:N_dot;
dotV2 = VD2*cos(1.5.*pi.*dot_vec2./N_dot)';
dotV2 = flipud(dotV2);
dot2 = spdiags(cat(1,zeros(N_tot-N_dot,1),dotV2),0,N_tot,N_tot);

NOT_dot = spdiags(cat(1,zeros(N_dot,1),ones(N_tot-2.*N_dot,1),zeros(N_dot,1)),0,N_tot,N_tot); % QD on both sides
%NOT_dot = spdiags(cat(1,zeros(N_dot,1),ones(N_tot-N_dot,1)),0,N_tot,N_tot); % QD on the left.
%NOT_dot = spdiags(cat(1,ones(N_tot-N_dot,1),zeros(N_dot,1)),0,N_tot,N_tot); % QD on the right.
%NOT_dot = spdiags(ones(N_tot,1),0,N_tot,N_tot); % No QD at all.

%%
H = zeros(4*N_tot,4*N_tot);
%SelfE = -lambda.*((omega)*t0s0 + Delta*txs0)./sqrt(Delta.^2 -(omega).^2); %%%%%
%H = kron(eyesm, (2*t - mu)*tzs0 + Vz*t0sx) + kron(barrier,Ebarrier*tzs0) + kron(dot1, tzs0) + kron(NOT_dot, SelfE) + kron(dot2,tzs0); % Diagonal Part --- Nanowire
H = kron(eyesm, (2*t - mu)*tzs0 + Vz*t0sx) + kron(barrier,Ebarrier*tzs0)+ kron(dot1, tzs0) + kron(NOT_dot, Delta*txs0) + kron(dot2,tzs0);

H = H + kron(band11sm, -t*tzs0) + kron(band1m1sm, 1i*alpha*tzsy); % Off-diagonal Part
end