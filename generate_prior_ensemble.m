function En=generate_prior_ensemble(n,J,Prior,Wall,Meas)
%Inputs: 
%%% n: n=log_2(N) where N is the number of cells used to discretise
%%% the computational domain of the wall's width J=ensemble size (i.e.
%%% number of ensemble members)
%%% Wall: structure with the geometry of the wall and time-stepping for the
%%% heat transfer model.
%%% Meas: structure that contains synthetic measurements of surface heat
%%% flux and near air temperatures, standard deviations of surface heat
%%% flux measurements, number of measurements within each assimilation subinterval
%%% Prior: structure that contains hyperparameters of the lognormal random
%%% fields for (k, c, T_0) and mean and standard deviation for (R_I, R_E)

%Outputs:
%%% En: structure with the prior ensemble (of size J) of (k, c, T_0, R_I, R_E).
%%% En.thermal{i} contains prior samples for the log(c) (i=1), log(k)
%%% (i=2), T_{0} (i=3), R_{I} (i=4) and R_{E} (i=5)
%%% En.thermal_means{i} and En.thermal_R{i} contains prior means and prior
%%% KL-based square root covariances.
%%% Samples from the prior are generated according to the algorithm
%%% described in Appendix A in De Simon, Iglesias, Jones, 
%%% and Wood (To appear in Energy and Buildings)



Nx=2^n;
Grid.Nx=Nx;
Grid.hx=Wall.L/Nx;
x_c=linspace(Grid.hx/2,Wall.L-Grid.hx/2,Grid.Nx)'; %define centers of cells
x_e=linspace(0,Wall.L,Grid.Nx+1)'; %define edges of cells

%Store prior means
En.thermal_means{1}=Prior.c_mean*ones(Nx,1); 
En.thermal_means{2}=Prior.k_mean*ones(Nx,1);
%Define prior mean for T_0 (see equation A.5 in Appendix)
T0_0=(Meas.T_int(1)-Meas.HF.int(1)*exp(Prior.R_I_mean+Prior.R_I_std^2/2));
T0_L=(Meas.T_ext(1)+Meas.HF.ext(1)*exp(Prior.R_E_mean+Prior.R_E_std^2/2));
En.thermal_means{3}=(T0_0+(T0_L-T0_0)*x_e/Wall.L);

En.thermal_means{4}=Prior.R_I_mean;
En.thermal_means{5}=Prior.R_E_mean;

%%% generate (squared root) covariance matrices for Gaussian random fields
En.thermal_R{1}=cova_matrix(x_c,Prior.c); 
En.thermal_R{2}=cova_matrix(x_c,Prior.k);
En.thermal_R{3}=cova_matrix(x_e,Prior.T0);

En.thermal_R{4}=Prior.R_I_std;
En.thermal_R{5}=Prior.R_E_std;
En.J=J;
En.n=n;
%%% generate the ensemble of Gaussian random fields
for en=1:En.J
    %%% prior ensemble of log(c)
    En.thermal{1}(:,en) =En.thermal_means{1}+En.thermal_R{1}*randn(Nx,1);
    %%% prior ensemble of log(k)
    En.thermal{2}(:,en)=En.thermal_means{2}+En.thermal_R{2}*randn(Nx,1);
    %%% prior ensemble of T_0
    En.thermal{3}(:,en)=En.thermal_means{3}+En.thermal_R{3}*randn(Nx+1,1);
end
%%% prior ensemble of R_I
En.thermal{4}=Prior.R_I_mean+Prior.R_I_std*randn(1,J);
%%% prior ensemble of R_E
En.thermal{5}=Prior.R_E_mean+Prior.R_E_std*randn(1,J);




