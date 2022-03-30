function rho=dosH_se_v14(t,Delta,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,V,s)
H = hse_v10(t,Delta,N_tot,alpha,mu,VD1,VD2,N_dot,Nbarrier,Ebarrier,Vz,lambda,V+1i*s);
temp = (V + 1i*s)*speye(4*N_tot);
G = inv(H - temp);
rho = (G - G')./(2i.*pi);