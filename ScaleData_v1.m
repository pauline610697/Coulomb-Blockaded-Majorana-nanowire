clear;
load Landauer_v55_ext_3.mat
%load Landauer_v56_L=150_Total_6_TwoQD_mu=2.5_lambda=1.4_VD1=1.0_VD2=4.0_Vc=4.2_skipBarrierSites_T=0.01.mat
%load Landauer_v55_L=150_Total_8_TwoQD_mu=2.5_lambda=1.4_VD1=1.0_VD2=4.0_Vc=4.2_skipBarrierSites_T=0.01.mat
Factor = 15/40;
N_even = 20;
N_odd = 21;
ngRange_even = (N_even - 1./2) - (N_even - ngRange - 1./2)./Factor; 
ngRange_odd = (N_odd - 1./2) - (N_odd - ngRange - 1./2)./Factor; 

[~,even_index_min] = min(abs(ngRange_even-19));
[~,even_index_max] = min(abs(ngRange_even-21));
[~,odd_index_min] = min(abs(ngRange_odd-19));
[~,odd_index_max] = min(abs(ngRange_odd-21));

G_Even_1 = G_Even(:,even_index_min:even_index_max);
G_Odd_1 = G_Odd(:,odd_index_min:odd_index_max);
ngStep_1 = ngStep./Factor;
ngRange_1 = 19:ngStep_1:21;

% Scale the conductance height: Make normal matal peak G=1, G in other regimes are relative to this.
% This is just for plotting purpose.
%Tn = (VzcNumber + VzNumber)/2; % Mean number of Vc and V_{z,max}
%pk_height = findpeaks(G_Even(Tn,:));
%modified_factor = 1/pk_height;
%modified_factor = 1/0.05325;
%modified_factor = 1/0.05;
%% One Copy
G_combine = G_Even_1(51:end,:) + G_Odd_1(51:end,:);
%G_combine = G_Even(101:end,:) + G_Odd(101:end,:);
%G_combine = G_combine.*modified_factor;
figure()
imagesc(VzRange(51:end),ngRange,G_combine')
colorbar
%caxis([0 2]) % Scale the color
title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.2$ meV, $\alpha=2.5$ meV, $\mu=2.5$ meV, $\lambda=1.4$ meV, $E_c=0.4$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=150$, $T=0.01$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$n_g$$','interpreter','latex','FontSize',16)
set(gca,'Ydir','normal')
%% Two Copy
ngRange_2 = 17:ngStep_1:21;
%ngRange_2 = 10:ngStep_1:50;
G_Even2 = cat(2,G_Even,G_Even(:,2:end));
G_Odd2 = cat(2,G_Odd,G_Odd(:,2:end));
%G_Even2 = cat(2,G_Even_1(101:end,:),G_Even_1(101:end,2:end));
%G_Odd2 = cat(2,G_Odd_1(101:end,:),G_Odd_1(101:end,2:end));
G_combine2 = G_Even2 + G_Odd2;
%G_combine2 = G_combine2.*modified_factor;
figure()
imagesc(VzRange,ngRange_2,G_combine2')
%imagesc(VzRange(101:end),ngRange_2,G_combine2')
colorbar
%caxis([0 0.1e-5]) % Scale the color
%xline(2.0,'--w');
%xline(2.6,'--w');
xline(sqrt(lambda^2+mu^2),'--w');
%xline(sqrt(Delta_0^2+mu^2),'--w');
xline(Vzc,'--w');
%caxis([0 0.07]) % Scale the color
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=2.5$ meV, $E_c=0.16$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.2$ meV, $\alpha=2.5$ meV, $\mu=2.5$ meV, $\lambda=1.4$ meV, $E_c=0.08$ meV, $V_{D1}=1.0$ meV, $l_D=26$, $L=150$, $T=0.01$ meV.','interpreter','latex','FontSize',16)
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.2$ meV, $\alpha=2.5$ meV, $\mu=2.5$ meV, $\lambda=1.4$ meV, $E_c=0.03$ meV, $V_{D1}=1.0$ meV, $V_{D2}=4.0$ meV, $l_D=26$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
xlabel('$\textsf{V}_\textsf{z}$  \sf{(meV)}','interpreter','latex','FontSize',20,'FontName','SansSerif')
ylabel('$\textsf{n}_\textsf{g}$','interpreter','latex','FontSize',20,'FontName','SansSerif')
%xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
%ylabel('$$n_g$$','interpreter','latex','FontSize',16)
set(gca,'Ydir','normal')
set(gca,'FontSize',20)
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',24)
%xticks(1:1:5)
yticks(17:0.5:21)
%% Plot Conductance
G_combine = G_Even + G_Odd;
%G_combine = G_combine.*modified_factor;
figure()
imagesc(VzRange,ngRange,G_combine')
colorbar
xline(sqrt(lambda^2+mu^2),'--w');
xline(Vzc,'--w');
%caxis([0 1]) % Scale the color
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=2.5$ meV, $E_c=3$ meV, $V_{D1}=2.0$ meV, $l_D=30$, $L=150$, $T=0.02$ meV.','interpreter','latex','FontSize',16)
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.2$ meV, $\alpha=2.5$ meV, $\mu=3.0$ meV, $\lambda=1.4$ meV, $E_c=3$ meV, $V_{D1}=0.7$ meV, $l_D=26$, $L=220$, $T=0.01$ meV.','interpreter','latex','FontSize',16)
%title('$t=25$ meV, $\Delta_0=0.9$ meV, $V_c=4.2$ meV, $\alpha=2.5$ meV, $\mu=3.0$ meV, $\lambda=1.4$ meV, $E_c=3$ meV, $V_{D1}=0.7$ meV, $V_{D2}=4.2$ meV, $l_D=26$, $L=150$, $T=0.01$ meV.','interpreter','latex','FontSize',16)
xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
ylabel('$$n_g$$','interpreter','latex','FontSize',16)
set(gca,'Ydir','normal')

%figure()
%imagesc(VzRange,ngRange_even,G_Odd_1')
%colorbar
%caxis([0 1]) % Scale the color
%title('$$G_{even}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=3.5$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=1.5$ meV, $E_c=3$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=80$, $T=0.02$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
%xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
%ylabel('$$n_g$$','interpreter','latex','FontSize',16)
%set(gca,'Ydir','normal')

%figure()
%imagesc(VzRange,ngRange_odd,G_Odd_1')
%colorbar
%caxis([0 1]) % Scale the color
%title('$$G_{odd}$$: $t=25$ meV, $\Delta_0=0.9$ meV, $V_c=3.5$ meV, $\alpha=2.5$ meV, $\mu=2.0$ meV, $\lambda=1.5$ meV, $E_c=0.17$ meV, $\rho_F=10^{-1}$, $V_{SC}=10^5$, $L=150$, $T=0.02$ meV, $s=10^{-3}$, threshold$=10^{-4}$.','interpreter','latex','FontSize',16)
%xlabel('$$V_z$$ (meV)','interpreter','latex','FontSize',16)
%ylabel('$$n_g$$','interpreter','latex','FontSize',16)
%set(gca,'Ydir','normal')
