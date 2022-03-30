%% Version 56-B-0: This code is to calculate dI/dV false plot/peak spacings as a function of ng and Vz for the multi-level system.
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
% Separate ls for electron & hole parts (use "LDOS_se_v3.m") ==> Difference from v14b
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
% Fix the prod(tanh(s)) issue by taking log(lambda) ==> Difference from v51
% Rewrite prod(tanh(s)) into lambda in Z to fix (1-lambda)=1 issue ==> Difference from v52
% Two QDs ==> Difference from v51-v53
% Soft gap modification, treat soft-gap region the same as above Vc ==> Difference from v55

clear;
tic;
C = 20; part = 0;
%% Parameters Setting
% Note that the length scale is in unit of lattice constant, which is 10nm.

load(['Landauer_v54_L=150_partA_',num2str(part),'.mat'])
% Basic Parameters
N = 20; % The minimum of N should be chosen as the # of energy levels we consider. i.e. N >= levelN.
vec_N = [(N-1) N (N+1)];

t = 25; %unit: meV
Delta_0 = 0.9; %unit: meV
Vzc = 4.2;
wireLength = 150; %unit: 10nm
alpha = 2.5; %unit: meV
mu = 2.5; %unit: meV
lambda = 1.4; %unit: meV
T = 0.01; %unit: meV
Ec = 3; %unit: meV

VD1 = 1; %unit: meV
VD2 = 4; %unit: meV
N_dot = 26; %unit: 10nm

Nbarrier = 2; %unit: 10nm
Ebarrier = 10; %unit: meV

N_tot = wireLength;
s = 1e-3; % Vstep = 1e-3; resolution is 5 times smaller than the width.
rho_F = 0.1;
V_SC = 10^5;
V_metal = 2.5;
threshold = 1e-2;

%% Construct the variable space
VzMin = 0; VzMax = 5.5; VzNumber = 551;
VzStep = (VzMax - VzMin)./(VzNumber - 1);
VzRange = linspace(VzMin,VzMax,VzNumber);

Vc2 = sqrt(lambda.^2 + mu.^2);
Vc2Number = floor(1 + (Vc2 - VzMin)./VzStep);

ngMin = 19; ngMax = 21; ngNumber = 2001; % (19.4; 19.6; 81)
ngStep = (ngMax - ngMin)./(ngNumber - 1);
ngRange = linspace(ngMin,ngMax,ngNumber);

Vmin = 0; Vmax = 1.1; Vnumber = 11001; % (-0.0025, 4, 1002)1101
Vstep = (Vmax - Vmin)./(Vnumber - 1); % = 1e-3
Vrange = linspace(Vmin,Vmax,Vnumber);

G_odd = zeros(VzNumber,ngNumber);
G_even = zeros(VzNumber,ngNumber);

G_odd_L = zeros(VzNumber,ngNumber);
G_odd_R = zeros(VzNumber,ngNumber);
G_even_L = zeros(VzNumber,ngNumber);
G_even_R = zeros(VzNumber,ngNumber);

level_num = zeros(1,C);
total_num = zeros(1,C);
E_less_T_num = zeros(1,C);
Q_0 = zeros(1,C);
Parity = zeros(1,C);

%%
parfor k = 1:C % C; VzNumber
    K = C.*part+k;
    Vz = VzMin + (K-1).*VzStep;
    disp(K);
    
    Delta1 = Delta_0.*sqrt(1 - (Vz./Vzc).^2).*(Vz<Vzc);
    %Delta1 = Delta_0;
    Delta2 = Delta_0.*sqrt(1 - (Vz./Vc2).^2).*(Vz<Vc2);
    
    if K<Vc2Number
        vecE = [firstE(k) dosmap1{k}];
    end
    
    Gamma_L = [];
    Gamma_R = [];
    Lambda_L = [];
    Lambda_R = [];
    vecE_temp = [];
    parity = 0;
    if K<Vc2Number
        locMin = locMin1{k};
        V_locMin = Vrange(locMin);
        Bound_N = bound_num(k); %%%
        for n = 1:Bound_N
           rho_w = 0;
           zero_div = 1;
           %disp(['bound state = ',num2str(n)]);

           if locMax1(k)<=50 && n==1 % zero-energy case
               zero_div = 2;
           end
           n1 = n;
           x1 = V_locMin(n1); % integration lower bound
           x2 = V_locMin(n1+1); % integration upper bound
                  
           x = x1:Vstep:x2; % integration resolution is Vstep.
           upperBound = (length(x)-1);
           for j = 1:upperBound
               %rho_w_1 = dosH_se_v14(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,x(j),s);
               %rho_w_2 = dosH_se_v14(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,x(j+1),s);               
               rho_w_1 = dosH_se_v16(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,x(j),s); % Soft Gap
               rho_w_2 = dosH_se_v16(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,x(j+1),s); % Soft Gap
               rho_w = rho_w + (rho_w_1 + rho_w_2);
           end
           rho_w = rho_w.*(x2-x1)./(2.*(length(x)-1)); % See the Trapezoidal Method from Matlab manual.
            
           [dege_WF,lambda_m] = eigs(rho_w,4.*N_tot,'LM');
           lambda_m = diag(lambda_m);
           lambda_m = real(lambda_m);
             
           index = find(lambda_m >= threshold); % 100.*
           %if n==1 && length(index)>1
           %    d = length(index)./zero_div;
           %else
           %    d = length(index);
           %end
           d = length(index);

           index = index(1:d);
           lambda_m = lambda_m(index); % selected d degeneracy
           dege_WF = dege_WF(:,index); % selected d degeneracy
           prefactor = sqrt(lambda_m'); % 1-by-d vector
           
           if n==1
               [Npsi_1,Nphi_1] = Majorana_rebuild(dege_WF(:,1),N_tot);
               H_BdG = kron(Npsi_1,Npsi_1') - kron(Nphi_1,Nphi_1');
               %Hami = hse_v10(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,0);
               Hami = hse_v11(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,0); % Soft gap
               Q_0(k) = FermionParityBdG(real(Hami + 10.*H_BdG));
               parity = (1-Q_0(k))./2;
               Parity(k) = parity;
               
               %if d==1
               %    dege_WF = Npsi_1;
               %else
               %    dege_WF = [Npsi_1 dege_WF(:,2:end)];
               %end
           end
                              
           dege_WF = prefactor.*dege_WF;                      
           vecE_temp = [vecE_temp repmat(vecE(n),1,d)];

           Gamma_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2);
           Gamma_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2);
           Lambda_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2);
           Lambda_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2);
            
           Gamma_L = [Gamma_L Gamma_l_temp];
           Gamma_R = [Gamma_R Gamma_r_temp];
           Lambda_L = [Lambda_L Lambda_l_temp];
           Lambda_R = [Lambda_R Lambda_r_temp];
        end
        
        tempE = Delta2 + 1e-16;
        metal1_E = [tempE];
        accuD = 0;
        n1 = n;
        while tempE < V_metal
            n1 = n1+1;
            %disp(['metal state = ',num2str(n1)]);
            
            %rho_w = dosH_se_v14(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,tempE,s);
            rho_w = dosH_se_v16(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,tempE,s); % Soft Gap
            [dege_WF,lambda_m] = eigs(rho_w,4.*N_tot,'LM');
            lambda_m = diag(lambda_m);
            lambda_m = real(lambda_m);
            
            index = find(lambda_m >= threshold);
            d = length(index);
            accuD = accuD + d;
            F = accuD./V_SC;
            if (tempE <= Delta1) && (tempE > Delta2)
                dEn = sqrt(Delta2.^2 + (accuD./(2.*rho_F.*V_SC)).^2) - tempE;
            elseif (tempE > Delta1)
                dEn = sqrt((Delta1.^2 + Delta2.^2 + (F./rho_F).^2).^2 - 4.*(Delta1.*Delta2).^2)./(2.*F./rho_F) - tempE;
            end
                
            lambda_m = lambda_m(index); % selected d degeneracy
            dege_WF = dege_WF(:,index); % selected d degeneracy
            prefactor = sqrt(lambda_m'.*dEn); % 1-by-(4.*N_tot) vecotr
            dege_WF = prefactor.*dege_WF;
            vecE_temp = [vecE_temp repmat(tempE,1,d)];
            tempE = tempE + dEn;
            metal1_E = [metal1_E tempE];
            
            Gamma_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2);
            Gamma_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2);
            Lambda_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2);
            Lambda_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2);
            
            Gamma_L = [Gamma_L Gamma_l_temp];
            Gamma_R = [Gamma_R Gamma_r_temp];
            Lambda_L = [Lambda_L Lambda_l_temp];
            Lambda_R = [Lambda_R Lambda_r_temp];
        end
        dosmap1{k} = [dosmap1{k} metal1_E(1:end-1)];
        level_num(k) = bound_num(k) + (length(metal1_E)-1);
    else 
        %DOS = arrayfun(@(V) dosH_se_v11(t,Delta1,N_tot,alpha,mu,VD,N_dot,Nbarrier,Ebarrier,Vz,lambda,V,s), Vrange);
        Q_0(k) = 1;
        tempE = Delta2 + 1e-16;
        metal2_E = [tempE];
        accuD = 0;
        n1 = 0;
        while tempE < V_metal
            n1 = n1+1;
            %disp(['metal state = ',num2str(n1)]);
            
            %rho_w = dosH_se_v14(t,Delta1,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,tempE,s);
            rho_w = dosH_se_v16(t,Delta1,Delta2,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,tempE,s); % Soft Gap
            [dege_WF,lambda_m] = eigs(rho_w,4.*N_tot,'LM');
            lambda_m = diag(lambda_m);
            lambda_m = real(lambda_m);
            
            index = find(lambda_m >= threshold);
            d = length(index);
            accuD = accuD + d;            
            F = accuD./V_SC;
            if (tempE <= Delta1) && (tempE > Delta2)
                dEn = sqrt(Delta2.^2 + (accuD./(2.*rho_F.*V_SC)).^2) - tempE;
            elseif (tempE > Delta1)
                dEn = sqrt((Delta1.^2 + Delta2.^2 + (F./rho_F).^2).^2 - 4.*(Delta1.*Delta2).^2)./(2.*F./rho_F) - tempE;
            end
            lambda_m = lambda_m(index); % selected d degeneracy
            dege_WF = dege_WF(:,index); % selected d degeneracy
            prefactor = sqrt(lambda_m'.*dEn); % (4.*N_tot)-by-1 vecotr
            dege_WF = prefactor.*dege_WF;
            vecE_temp = [vecE_temp repmat(tempE,1,d)];
            tempE = tempE + dEn;
            metal2_E = [metal2_E tempE];
            
            Gamma2_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2);
            Gamma2_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2);
            Lambda2_l_temp = not(parity).*(abs(dege_WF(4.*Nbarrier+1,:)).^2 + abs(dege_WF(4.*Nbarrier+2,:)).^2) + parity.*(abs(dege_WF(4.*Nbarrier+3,:)).^2 + abs(dege_WF(4.*Nbarrier+4,:)).^2);
            Lambda2_r_temp = not(parity).*(abs(dege_WF(4*N_tot - 3,:)).^2 + abs(dege_WF(4*N_tot - 2,:)).^2) + parity.*(abs(dege_WF(4*N_tot - 1,:)).^2 + abs(dege_WF(4*N_tot,:)).^2);
            
            Gamma_L = [Gamma_L Gamma2_l_temp];
            Gamma_R = [Gamma_R Gamma2_r_temp];
            Lambda_L = [Lambda_L Lambda2_l_temp];
            Lambda_R = [Lambda_R Lambda2_r_temp];
        end 
        level_num(k) = length(metal2_E)-1;
    end
    vecE = vecE_temp;
    temp = tanh(vecE./(2.*T));    
    
    total_num(k) = length(vecE);
    Total_N = total_num(k);
    y = zeros(1,Total_N);
    
    x1 = vecE(1)./(2.*T);
    y(1) = 2.*x1 + log((1+exp(-2.*x1))./2);    
    for p = 2:Total_N
        x = vecE(p)./(2.*T);
        z = 2.*x + log((1+exp(-2.*x))./2);
        a = max([y(p-1) z]);
        b = min([y(p-1) z]);
        y(p) = b - log(1 + exp(-(a-b)) - exp(-a));
    end
    gamma = exp(-y(Total_N));
    prod_T = 1-gamma;
    prod_T_noP = prod_T./temp;
    
    temp_index = find(vecE<T);
    E_less_T_num(k) = length(temp_index);
    
    for m = 1:ngNumber
        ng = ngMin + (m-1).*ngStep;
        U_N = Ec.*(vec_N-ng).^2;
        dU = diff(U_N); % dU = U(N) - U(N-1)
        
        if mod(N,2)==0 % Make dU(1) is for N odd and dU(2) is for N even.
            dU = fliplr(dU);
        end 
        
        Z_odd = (1 + exp(-vecE./T)).*(gamma + exp(-(-Q_0(k).*dU(1)./T)).*(2-gamma));
        Z_even = (1 + exp(-vecE./T)).*(gamma + exp(-(Q_0(k).*dU(2)./T)).*(2-gamma));
        tempM_odd_L = zeros(1,total_num(k));
        tempM_odd_R = zeros(1,total_num(k));
        tempM_even_L = zeros(1,total_num(k));
        tempM_even_R = zeros(1,total_num(k));
        for n = 0:1
            for Q = -1:2:1 % Q is either (-1) or (+1)
                dE_odd = vecE.*(1 - 2.*n) + Q.*dU(1); % N is odd.
                dE_even = vecE.*(1 - 2.*n) - Q.*dU(2); % N is even.
                
                Fp_odd = 0.5.*exp(-(vecE - Q_0(k).*dU(1))./(2.*T)).*sech(dE_odd./(2.*T))./T.*(1 + Q.*Q_0(k).*((-1).^n).*prod_T_noP)./Z_odd;
                Fp_even = 0.5.*exp(-(vecE + Q_0(k).*dU(2))./(2.*T)).*sech(dE_even./(2.*T))./T.*(1 + Q.*Q_0(k).*((-1).^n).*prod_T_noP)./Z_even;
                
                tempM_odd_L = tempM_odd_L + Fp_odd.*((1-n).*Gamma_L + n.*Lambda_L);
                tempM_odd_R = tempM_odd_R + Fp_odd.*((1-n).*Gamma_R + n.*Lambda_R);
                tempM_even_L = tempM_even_L + Fp_even.*((1-n).*Gamma_L + n.*Lambda_L);
                tempM_even_R = tempM_even_R + Fp_even.*((1-n).*Gamma_R + n.*Lambda_R);
            end
        end
        
        G_odd_L(k,m) = sum(tempM_odd_L); % Sum over the changed energy states and configurations
        G_odd_R(k,m) = sum(tempM_odd_R); % Sum over the changed energy states and configurations
        G_even_L(k,m) = sum(tempM_even_L); % Sum over the changed energy states and configurations
        G_even_R(k,m) = sum(tempM_even_R); % Sum over the changed energy states and configurations
        
        G_odd(k,m) = G_odd_L(k,m).*G_odd_R(k,m)./(G_odd_L(k,m) + G_odd_R(k,m));
        G_even(k,m) = G_even_L(k,m).*G_even_R(k,m)./(G_even_L(k,m) + G_even_R(k,m));        
    
        disp(['(Vz,ng)=(',num2str(Vz),',',num2str(ng),')']);
    end
end
toc;
save(['Landauer_v56_L=150_partB_',num2str(part),'.mat'])