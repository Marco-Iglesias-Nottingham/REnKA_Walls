function [Meas, Truth]=generate_synthetic_data(noise,Wall)
%Inputs: 
%%%% noise: standard deviation of the noise added to the noise-free
%%% measurements of heat flux.
%%% Wall: structure with the geometry of the wall and
%%% time-stepping for the heat transfer model.
%Outputs:
%%% Meas: structure that contains synthetic measurements of surface heat 
%%% flux and near air temperatures, standard deviations of surface heat 
%%% flux measurements, number of measurements within each assimilation subinterval
%%% Truth: structure that contains the true thermophysical
%%% properties used to generate synthetic data (see section 3.1 from
%%% ``Quantifying uncertainty in thermophysical properties of walls 
%%% by means of Bayesian inversion'' by De Simon, Iglesias, Jones, 
%%% and Wood (To appear in Energy and Buildings)


n=9;
Nx=2^n; %number of discretisation points
hx=Wall.L/Nx;


%Define true synthetic measurements of temperature:
rng('default')
rng(12689);


%%%%%%%%%%%%%%%% Generate stochastic (Gaussian) process that we use to 
%%perturb the internal and external temperatures
%%%%%%%%%%%%%%%% temperatures
tt=linspace(0,3*Wall.T,3*Wall.Nt);
T(1)=sqrt(2); T(2)=1.25; T(3)=Wall.T/30;
R_t=cova_matrix(tt,T);
t_int=R_t*randn(3*Wall.Nt,1);
t_ext=R_t*randn(3*Wall.Nt,1);

%%% Define the internal and external synthetic measurements defined on an
%%% interval of length [-T,2T] 
%%% recall we use [-T,0] to generate initial temperature distribution
%%% Assimilation is done in [0,T] and predictions are made for [T, 2T]

T_int=(20+1.0*sin(tt*Wall.T*2))+t_int';
T_ext=(8.5+2.5*cos(tt*Wall.T))+t_ext';



%Define true thermophysical properties:
%x=linspace(0,Wall.L,Nx)';
x=linspace(hx/2,Wall.L-hx/2,Nx)';  

%%True specific heat capacitance 
c_truth=log(10^6*0.7)*ones(Nx,1);
Truth.thermal{1}=c_truth+0.40.*((-0.05<x)&(x<0.0443*0.7))+1.35.*((0.0443*0.7<x)&(x<0.0443*1.0))...
    +0.4.*((0.0443*1.0<x)&(x<0.0443*1.9))-1.35*((0.0443*1.9<x)&(x<0.0443*2.4))+0.40.*((0.0443*2.4<x)&(x<0.0443*3.8))...'
    +1.35*((0.0443*3.8<x)&(x<0.0443*4.4))+0.4*((0.0443*4.4<x)&(x<0.0443*5.8)) -1.5.*((0.0443*5.8<x)&(x<0.0443*6.2))...'
    +0.40.*((0.0443*6.2<x)&(x<0.0443*8.6));

Truth.exp_thermal{1}=exp(Truth.thermal{1});

%%True thermal conductivity
u_truth=log(0.5)*ones(Nx,1);%
Truth.thermal{2}=u_truth+0.95.*((-0.05<x)&(x<0.0443*0.7))+1.5.*((0.0443*0.7<x)&(x<0.0443*1.0))...
    +0.95.*((0.0443*1.0<x)&(x<0.0443*1.9))-1.0*((0.0443*1.9<x)&(x<0.0443*2.4))+0.95.*((0.0443*2.4<x)&(x<0.0443*3.8))...'
    +1.5*((0.0443*3.8<x)&(x<0.0443*4.4))+0.95*((0.0443*4.4<x)&(x<0.0443*5.8))- 1.0.*((0.0443*5.8<x)&(x<0.0443*6.2))...'
    +0.95.*((0.0443*6.2<x)&(x<0.0443*8.6));

Truth.exp_thermal{2}=exp(Truth.thermal{2});

%True internal surface resistance:
Truth.thermal{4}=0.13;
%True external surface resistance
Truth.thermal{5}=0.04;
%%True C-value:
Truth.thermal{6}=sum(exp(Truth.thermal{1})*hx);
%%True U-value:
Truth.thermal{7}=1/(sum(1./Truth.exp_thermal{2})*hx+Truth.thermal{4}+Truth.thermal{5});

%%True initial temperature distribution
x_t=linspace(0,Wall.L,Nx+1);
%%start with a linear profile that interpolates T_int and T_ext at time
%%t=-T
T0=T_int(1)-(T_int(1)-T_ext(1))*x_t/Wall.L;
%%then use this T0 to run the Heat Equation for an interval of time [-T,0]
Truth.thermal{3}=Heat_FEM(hx,Nx,Wall.Nt,Wall,exp(Truth.thermal{2}),...'
    exp(Truth.thermal{1}),T0',Truth.thermal{4},Truth.thermal{5},T_int(1:Wall.Nt),T_ext(1:Wall.Nt),0);
%%Define the initial true temperature as the final temperature at t=0:





Meas.T_all.int=T_int;  %save all measurements of near air temperatures in the interval [-T,2T]
Meas.T_all.ext=T_ext;
Meas.T_int=T_int(Wall.Nt+1:end); %save measurements of near air temperatures in the interval [0,2T]
Meas.T_ext=T_ext(Wall.Nt+1:end);
Meas.DeltaT=30;  %number of observations within the assimilation subinterval [t_{m-1},t_{m}]

%%Generate synthetic heat fluxes on the interval [0,2T]
HF=Heat_FEM(hx,Nx,2*Wall.Nt,Wall,exp(Truth.thermal{2}),exp(Truth.thermal{1}),...'
    Truth.thermal{3},Truth.thermal{4},Truth.thermal{5},Meas.T_int,Meas.T_ext,1);
%%generate noise
gauss_noise_int=randn(2*Wall.Nt,1);
gauss_noise_ext=randn(2*Wall.Nt,1);


MM=floor(2*Wall.Nt/Meas.DeltaT);
for m=1:MM
    ind_m_1=1+(m-1)*Meas.DeltaT;
    ind_m=m*Meas.DeltaT;
    %%% compute measurement standard deviation on each subinterval[t_{m-1},t_{m}]
    %%% see equation (25) from section 3.1 from De Simon et-al
    Meas.stda.int(ind_m_1:ind_m)=noise*sum(abs(HF.int(ind_m_1:ind_m)))/Meas.DeltaT;
    Meas.stda.ext(ind_m_1:ind_m)=noise*sum(abs(HF.ext(ind_m_1:ind_m)))/Meas.DeltaT;
end

%%% synthetic measurements of surface heat flux (contaminated with Gaussian noise)
Meas.HF.int=HF.int+Meas.stda.int'.*gauss_noise_int; 
Meas.HF.ext=HF.ext+Meas.stda.ext'.*gauss_noise_ext;



