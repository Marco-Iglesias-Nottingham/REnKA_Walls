clear all
close all
%%%Main file to reproduce some of the results from "Quantifying uncertainty
%%in thermophysical properties of walls by means of Bayesian inversion" by 
%%De Simon, Iglesias, Jones and Wood (To appear in Energy and Buildings). 


%Define wall parameters and time-stepping for heat transfer model
Wall.L=0.31;  %Walls width [0,Wall.L]
%frequency of observations (also used as time stepping for the 1D heat
%equation
Wall.tau=60*5; %5 mins time step
Wall.Nt=1800;% total number of observation within the assimilation interval.
Wall.T=Wall.Nt*Wall.tau/60/60/24;%length of assimilation interval 6.25 days

noise=0.05; %%% define noise level in the heat flux observation (e.g. 5%)


%%% Generate synthetic measurements (surface heat fluxes and near air
%%% temperatures)

[Meas, Truth]=generate_synthetic_data(noise,Wall);

%%% Note: 
%%% If real data is used comment the above line and (i) include internal and 
%%% external near air temperatures in Meas.T_int, Meas.T_ext, (ii) include
%%% surface heat flux measurements in Meas.HF.int and Meas.HF.ext and (iii)
%%% compute heat flux measurement standard deviations and include them in 
%%% Meas.stda.int and Meas.stda.ext (e.g. via equation (25)). The code
%%% assumes all these measurements are collected at the same
%%% observation/monitoring times at a frequency specified by Wall.tau. 
%%% Include in Meas.DeltaT;  the number of observations that you wish to
%%% assimilate on each assimilation subinterval [t_{m-1},t_{m}]


J=1e3; %% Ensemble size (number of particles)

%%% Define Prior ensemble parameters

%%% lognormal Gaussian random fields for k, c and T_0
%%% See Appendix A (equation A.3) from De Simon et-al)
%%% The parameters below are required to specify Whittle-Mattern
%%% correlation function

Prior.k(1)=0.65;  %sigma
Prior.k(2)=1.05;  %nu  
Prior.k(3)=Wall.L/50; %l 
Prior.k_mean=log(0.75); %omega (see equation A.2)

Prior.c(1)=0.7;
Prior.c(2)=1.05;  
Prior.c(3)=Wall.L/50; 
Prior.c_mean=log(7.5e5);

Prior.T0(1)=sqrt(3.5); 
Prior.T0(2)=1.05; 
Prior.T0(3)=Wall.L/30;
%mean (i.e. omega) for T_0 is defined within generate_prior_ensemble below)


%Define mean and standard deviation for the log-normal priors for R_I and
%R_E

Prior.R_I_mean=log(0.1);
Prior.R_I_std=0.5;
Prior.R_E_mean=log(0.07);
Prior.R_E_std=0.5;


n=7; %%% 2^n defines the spatial discretisation of the wall's width [0,Wall.L].
%%% In order to avoid "inverse crimes" here we use different spatial 
%%% discretisation here than the one used for the generation of synthetic data (n=9). 

%%% generate initial ensemble 
En=generate_prior_ensemble(n,J,Prior,Wall,Meas);

%%% Run REnKA to compute approximate posteriors of input parameters 
%%% k, c, T_0, R_{int}, R_{ext} and to compute predictive distributions of
%%% heat fluxes on the predictive interval [T,2T]
Pos=REnKA(En,Wall,Meas);

T=floor(Wall.Nt/Meas.DeltaT); %Final assimilation time

%%% visualise synthetic data and statistics from the posterior.
visualise_results(Meas,Wall,Pos,En,Truth,T)


