%% Version 54-A-0: This code is to calculate dI/dV false plot/peak spacings as a function of ng and Vz for the multi-level system.
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
