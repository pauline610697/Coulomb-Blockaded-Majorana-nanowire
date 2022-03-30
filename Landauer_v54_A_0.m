%% Version 54-A-0: This code is to calculate dI/dV false plot/peak spacings as a function of ng and Vz for the multi-level system.
%  Using Landauer Formula. Phys. Rev. Lett. 68.2512, Eq.(13)
%  Superconducting gap change with the magnetic field
%  Include self-energy
%  Vectorization
%  Use Jay's simplification: only need to sum over changed orbitals p.
%  Deal with the degeneracy issue. => Difference between v11 and v12.
%  Find the root of the determinant(H_nw + selfE - \omega)==0.
% => Count the number of negative eigenstates, and find the difference \omega points.
% Include metallic part energy levels.
% Use LDOS to find WF
% Normalize the LDOS by using "LDOS_se_v2.m" ==> Difference from v14a
% Separate DOS for electron & hole parts (use "LDOS_se_v3.m") ==> Difference from v14b
% Include states above SC gap ==> Difference from v14c
% Use SC bulk gap DOS above SC gap ==> Difference from v14d
% eigen-energy is modified with 1st-order perturbation from the self-energy ==> Difference from v21.
% Deal with the degenerate wave-function ==> Difference from v22. (See 20191227 Jay's notes.)
% Use the exact expression for the density matrix in Jay's notes ==> Use dosH_se_v7.m ==> Difference from v25.
% Only count "distinct" energy levels to do the degeneracy. ==> Difference from v26.
% Find the eigenvalues with threshold 0.01. ==> Difference from v27
% Combination of v28 and v20. ==> Difference from v28
% Group and ensemble average the DOS peaks ==> Difference from v29
%  Use simplified F factor (cancellation with Z) ==> Difference from v31
% Do not fix the number of energy levels to be V_metal_N ==> Difference from v32
% Count the number of states with energy less than the temperature. ==> Difference from v33
% 1. rho_w -> rho_w(k,n) ==> Difference from v34
% 2. above the gap: density matrix multiplying (delta E_n) firstly ==> Difference from v34
% 3. Fix the n-dependent degneracy issue ==> Difference from v34
% Get rid of first zero energy; get rid of the case WITHOUT bound state ==> Difference from v38
% Different approach for picking the energy level above the gap ==> Difference from v39
% Same method as v40, but with QD
% Not saving so many variables ==> Difference from v42
% Add Jay's Pfaffian factor for zero-crossing ==> Difference from v43
% Interchange n with (1-n), Q with -Q, and Gamma with Lambda. ==> Difference from v45
% Consider MORE THAN the first energy level ==> Difference from v49
% Add parity when choosing the wave-function
% Rebuild zero-energy degenerate states from Majorana basis ==> Difference from v48 and v48b
% Unify F_p(n,Q) ==> Difference from v50
% Two QDs ==> Difference from v51-v53

clear;
tic;
C = 20; part = 0;
%% Parameters Setting
% Note that the length scale is in unit of lattice constant, which is 10nm.

% Basic Parameters
N = 20; % The minimum of N should be chosen as the # of energy levels we consider. i.e. N >= levelN.

t = 25; %unit: meV
Delta_0 = 0.9; %unit: meV
Vzc = 4.2;
wireLength = 150; %unit: 10nm
alpha = 2.5; %unit: meV
mu = 2.5; %unit: meV
lambda = 1.4; %unit: meV
T = 0.01; %unit: meV
Ec = 3; %unit: meV

VD1 = 4; %unit: meV
VD2 = 4; %unit: meV
N_dot = 26; %unit: 10nm

Nbarrier = 0; %unit: 10nm
Ebarrier = 0; %unit: meV

N_tot = wireLength;
s = 1e-3; % Vstep = 1e-3; resolution is 5 times smaller than the width.
rho_F = 0.1;
V_SC = 10^5;
V_metal = 2.5;
threshold = 1e-4;
%% Construct the Hamiltonian
VzMin = 0; VzMax = 5.5; VzNumber = 551;
VzStep = (VzMax - VzMin)./(VzNumber - 1);
VzRange = linspace(VzMin,VzMax,VzNumber);

VzcNumber = round(1 + (Vzc - VzMin)./VzStep);

Ng_even = zeros(1,VzNumber);
Ng_odd = zeros(1,VzNumber);

Vmin = 0; Vmax = 1.1; Vnumber = 11001; % (-0.0025, 4, 1002)1101
Vstep = (Vmax - Vmin)./(Vnumber - 1); % = 1e-3
Vrange = linspace(Vmin,Vmax,Vnumber);

DOS = zeros(1,Vnumber);
%DOS = zeros(VzNumber,Vnumber);
firstE = zeros(1,C); % (VzcNumber-1)
dosmap1 = cell(1,C);
bound_num = zeros(1,C);
locMin1 = cell(1,C); % index vector of local minimums for each Vz
locMax1 = zeros(1,C);

%%
parfor k = 1:C %%%%% (VzcNumber-1)
    K = C.*part+k;
    Vz = VzMin + (K-1).*VzStep;
    disp(K);
    
    Delta1 = Delta_0.*sqrt(1 - (Vz./Vzc).^2).*(Vz<Vzc);
    %Delta1 = Delta_0;
    %Vc2 = sqrt(lambda.^2 + mu.^2);
    %Delta2 = Delta_0.*sqrt(1 - (Vz./Vc2).^2).*(Vz<Vc2);
    
    DOS = arrayfun(@(V) dosH_se_v13(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,V,s), Vrange);
    %DOS = arrayfun(@(V) dosH_se_v15(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,V,s), Vrange);
    [~,locMax] = findpeaks([0 DOS]);
    locMax1(k) = locMax(1); % Save out
    [~,locMin] = findpeaks(-DOS);
    w_peak = [];
   
    locMin = [1 locMin];
    bound_num(k) = length(locMin)-1;
    locMin1{k} = locMin; % Save out    
    V_locMin = Vrange(locMin);
    rho = DOS;
    
    for n = 1:bound_num(k)
         x1 = V_locMin(n); index_1 = locMin(n);
         x2 = V_locMin(n+1); index_2 = locMin(n+1);
         x = x1:Vstep:x2;
         rho_n = rho(index_1:index_2);
         upper = trapz(x,rho_n.*x);
         lower = trapz(x,rho_n);
         AvgE = upper./lower;
         if AvgE<1.05*Delta1 % Change Delta1 to Delta2 for Soft-gap.
             w_peak = [w_peak AvgE];
         else
             bound_num(k) = bound_num(k)-1;
         end
    end
    
    if length(w_peak)>1
        firstE(k) = w_peak(1);
        dosmap1{k} = w_peak(2:end);
    elseif length(w_peak)==1
        firstE(k) = w_peak(1);
    else
        firstE(k) = 1e-16;
    end      
end
save(['Landauer_v54_L=150_partA_',num2str(part),'.mat'])
toc;