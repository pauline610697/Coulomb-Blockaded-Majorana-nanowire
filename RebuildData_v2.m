% Rebuild data for "2"-- no QD, with self-energy, Zeeman-field-varying SC gap
clear;

N = 20;
VzNumber = 551; ngNumber = 2001;
G_Even = zeros(VzNumber,ngNumber);
G_Odd = zeros(VzNumber,ngNumber);
Ng_Even = zeros(1,VzNumber);
Ng_Odd = zeros(1,VzNumber);

C = 20; % 20
%%
for i = 0:26 % 26
    %name = sprintf(['Landauer_v50_',num2str(i),'_partB.mat']);
    name = sprintf(['Landauer_v55_L=150_partB_',num2str(i),'.mat']);
    load(name);
    G_Even(C.*i+1:C.*(i+1),:) = G_even(1:C,:);
    G_Odd(C.*i+1:C.*(i+1),:) = G_odd(1:C,:);
    
    [MaxValue1,index1] = max(G_even(1:C,:),[],2);
    Ng_Even(C.*i+1:C.*(i+1)) = ngMin + (index1 - 1).*ngStep;
    
    [MaxValue2,index2] = max(G_odd(1:C,:),[],2);
    Ng_Odd(C.*i+1:C.*(i+1)) = ngMin + (index2 - 1).*ngStep;
end
i = 27; r1 = 11; % 27; 11
name = sprintf(['Landauer_v55_L=150_partB_',num2str(i),'.mat']);
load(name);
G_Even(C.*i+1:C.*i+r1,:) = G_even(1:r1,:);
G_Odd(C.*i+1:C.*i+r1,:) = G_odd(1:r1,:);
    
[MaxValue1,index1] = max(G_even(1:r1,:),[],2);
Ng_Even(C.*i+1:C.*i+r1) = ngMin + (index1 - 1).*ngStep;
    
[MaxValue2,index2] = max(G_odd(1:r1,:),[],2);
Ng_Odd(C.*i+1:C.*i+r1) = ngMin + (index2 - 1).*ngStep;

%for i = 24:26 % 26
    %name = sprintf(['Landauer_v50_',num2str(i),'_partB.mat']);
%    name = sprintf(['Landauer_v55_L=150_partB_',num2str(i),'.mat']);
%    load(name);
%    G_Even(C.*i+1-r1:C.*(i+1)-r1,:) = G_even(1:C,:);
%    G_Odd(C.*i+1-r1:C.*(i+1)-r1,:) = G_odd(1:C,:);
    
%    [MaxValue1,index1] = max(G_even(1:C,:),[],2);
%    Ng_Even(C.*i+1-r1:C.*(i+1)-r1) = ngMin + (index1 - 1).*ngStep;
    
%    [MaxValue2,index2] = max(G_odd(1:C,:),[],2);
%    Ng_Odd(C.*i+1-r1:C.*(i+1)-r1) = ngMin + (index2 - 1).*ngStep;
%end

%i = 27; r2 = 11; % 27; 11
%name = sprintf(['Landauer_v55_L=150_partB_',num2str(i),'.mat']);
%load(name);
%G_Even(C.*i+1-r1:C.*i+r2-r1,:) = G_even(1:r2,:);
%G_Odd(C.*i+1-r1:C.*i+r2-r1,:) = G_odd(1:r2,:);

%[MaxValue1,index1] = max(G_even(1:r2,:),[],2);
%Ng_Even(C.*i+1-r1:C.*i+r2-r1) = ngMin + (index1 - 1).*ngStep;
    
%[MaxValue2,index2] = max(G_odd(1:r2,:),[],2);
%Ng_Odd(C.*i+1-r1:C.*i+r2-r1) = ngMin + (index2 - 1).*ngStep;
    
if mod(N,2)==0
    Se = Ng_Odd - Ng_Even;
    So = (Ng_Even + 2) - Ng_Odd;
else
    Se = (Ng_Odd + 2) - Ng_Even;
    So = Ng_Even - Ng_Odd;
end

save Landauer_v55_ext_3--.mat
%save Landauer_v55_L=150_Total_8_TwoQD_lambda=2.5_VD1=2.0_skipBarrierSites.mat
%% Plot Conductance
figure()
imagesc(VzRange,ngRange,G_Even')
colorbar
%caxis([0 1]) % Scale the color
%title('$$G_{even}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=1.5$ meV, $E_c=3$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=80$, $T=0.02$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
title('$$G_{even}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=4.0$ meV, $\lambda=2.5$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%title('$$G_{even}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$n_g$$','interpreter','latex','FontSize',16)
set(gca,'Ydir','normal')

figure()
imagesc(VzRange,ngRange,G_Odd')
colorbar
%caxis([0 1]) % Scale the color
%title('$$G_{odd}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=1.5$ meV, $E_c=3$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=80$, $T=0.02$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
title('$$G_{odd}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=4.0$ meV, $\lambda=2.5$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%title('$$G_{odd}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$n_g$$','interpreter','latex','FontSize',16)
set(gca,'Ydir','normal')

%% Plot Conductance peak spacings
figure()
plot(VzRange,Se,'.','DisplayName','N is even.');
%plot(VzRange(1:470),Se(1:470),'.','DisplayName','N is even.');
hold on;
%figure()
plot(VzRange,So,'.','DisplayName','N is odd.');
%plot(VzRange(1:470),So(1:470),'.','DisplayName','N is odd.');

%title('$t=25$ meV, $\Delta_0=0.9$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=1.5$ meV, $E_c=3$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=80$, $T=0.02$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=4.0$ meV, $\lambda=2.5$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.7$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%axis([VzMin VzMax -1.5 1.5])
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$n_g$$','interpreter','latex','FontSize',16)