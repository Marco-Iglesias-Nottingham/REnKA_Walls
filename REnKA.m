function Pos=REnKA(En,Wall,Meas)

%Inputs: 
%%% En: structure with the prior ensemble (of size J) of (k, c, T_0, R_I, R_E).
%%% Wall: structure with the geometry of the wall and time-stepping for the
%%% heat transfer model.
%%% Meas: structure that contains synthetic measurements of surface heat
%%% flux and near air temperatures, standard deviations of surface heat
%%% flux measurements, number of measurements within each assimilation subinterval

%Outputs:
%%% Pos: structure with some statistics of the posterior ensemble of (k, c,
%%% T_0, R_I, R_E) computed via the ensemble Kalman inversion algorithm
%%% described in Appendix B of "Quantifying uncertainty in thermophysical 
%%% properties of walls by means of Bayesian inversion" by De Simon, 
%%% Iglesias, Jones and Wood (To appear in Energy and Buildings). 
%%% Pos includes posterior credible intervals and posterior mean of
%%% (k, c, T_0) computed at each assimilation time t_{m}. Pos also stores 
%%% the full ensembles of R_I, R_E, C-value and U-value. It also
%%% computes credible intervals and mean of the predictive distributions
%%% after assimilation (i.e. in the predictive interval [T,2T]. In
%%% addition, Pos provides chi-square and average interval score of
%%% internal and external surface heat fluxes over the predictive interval 


Nx=2^En.n;

Grid.Nx=Nx;
Grid.hx=Wall.L/Nx;

J=En.J;
OneN=ones(J);
OneN=(1/J)*OneN;

%Pos is a structure that will store all posterior quantities of interest
Pos.Predictive.mean=[];
Pos.Predictive.low=[];
Pos.Predictive.high=[];



M_final=floor(Wall.Nt/Meas.DeltaT); %Final assimilation time

J_thresh=J/3; %parameter for REnKA (See B.5 and Algorithm 2)

for m=1:M_final
    clc
    ind_m_1=1+(m-1)*Meas.DeltaT;
    ind_m=m*Meas.DeltaT;
    disp(['Assimilating measurements between: ', ...'
        num2str(ind_m_1*Wall.tau/60/60/24,'% 10.2f'), '  and  ', num2str(ind_m*Wall.tau/60/60/24,'% 10.2f'), '  days'])
    pause(0.1)
    %%% collect measurements in the interval [t_{m-1},t_{m}]
    meas_m=[Meas.HF.int(ind_m_1:ind_m);Meas.HF.ext(ind_m_1:ind_m)];
    %%% Define measurement error covariance matrices
    Gamma=diag([Meas.stda.int(ind_m_1:ind_m),Meas.stda.ext(ind_m_1:ind_m)].^2);
    inv_sqrt_Gamma=diag([1./(Meas.stda.int(ind_m_1:ind_m)),1./(Meas.stda.ext(ind_m_1:ind_m))]);
    sqrt_Gamma=diag([Meas.stda.int(ind_m_1:ind_m),Meas.stda.ext(ind_m_1:ind_m)]);
    
    for i=1:5
        if (i~=3)
            En.exp_thermal{i}=exp(En.thermal{i});
        end
    end
    En.exp_thermal{3}=En.thermal{3};
    %%% Compute and store percentiles for (k,c,T_0)
    for i=1:3
        Pos.thermal{i}.stat{1}(:,m) = prctile(En.exp_thermal{i}',2.5);
        Pos.thermal{i}.stat{2}(:,m) = prctile(En.exp_thermal{i}',97.5);
        Pos.thermal{i}.stat{3}(:,m) = mean(En.exp_thermal{i}');
    end
    %%Store the ensemble of R_I, R_E, C-value and U_value
    Pos.thermal{4}(:,m)=En.exp_thermal{4};
    Pos.thermal{5}(:,m)=En.exp_thermal{5};
    Pos.thermal{6}(:,m)=Grid.hx*sum(En.exp_thermal{1})';
    Pos.thermal{7}(:,m)=(1./(Grid.hx*sum(1./En.exp_thermal{2})+En.exp_thermal{4}+En.exp_thermal{5}))';
    Pos.En=En;
    %%% set parameters for REnKA (See Appendix B, Algorithm 3)
    iter=0;
    phi_perm_r=0;
    E=[];
    D=[];
    while (phi_perm_r<1)
        iter=iter+1;
        phi_perm_r_1=phi_perm_r;
        for en=1:J
            %%% run heat flux simulator from time 0 to time t_m:
            HF=Heat_FEM(Grid.hx,Grid.Nx,m*Meas.DeltaT,Wall,exp(En.thermal{2}(:,en)),...'
                exp(En.thermal{1}(:,en)),En.thermal{3}(:,en),...'
                exp(En.thermal{4}(:,en)),exp(En.thermal{5}(:,en)),Meas.T_int,Meas.T_ext,1);
            Heat_m(:,en)=[HF.int(ind_m_1:ind_m);HF.ext(ind_m_1:ind_m)]; %get fluxes on interval [t_{m-1},t_{m}]
            Z(:,en)=inv_sqrt_Gamma*(Heat_m(:,en) - meas_m); %compute weighted data misfit
            log_W(en)=1/2*norm(Z(:,en))^2; %compute weights
        end
        phi_perm_r=find_temp(J_thresh,log_W,phi_perm_r_1); % get tempering parameters
        alpha=1/(phi_perm_r-phi_perm_r_1);  %REnKA regularisation parameter
        for jj=1:J
            E(:,jj)=sqrt(alpha)*sqrt_Gamma*randn(length(meas_m),1); %noise perturbation
        end
        meanE=E*OneN;
        E=E-meanE;
        for jj=1:J
            D(:,jj)=meas_m + E(:,jj); %add noise to measurements
        end
        %%% compute covariances and cross-covariances for REnKA:
        mean_flux=Heat_m*OneN;
        D_flux=Heat_m-mean_flux;
        C_flux_flux=1/(J-1)* (D_flux*D_flux');
        omega=(C_flux_flux +alpha*Gamma)\(D-Heat_m);
        %%% Update parameters (log(c),log(k),T_0,log(R_I),log(R_E))
        for i=1:5
            En.thermal{i}=En.thermal{i}+ 1/(J-1)*(En.thermal{i}-En.thermal_means{i})*D_flux'*omega;
            En.thermal_means{i}=En.thermal{i}*OneN;
        end
        
    end
    %%% end of REnKA with tempering parameter phi_r=1
    
end

%%% Recompute empirical covariance matrices and SVD for resampling 
for i=1:5
    En.thermal_means{i}=En.thermal{i}*OneN;
    En.thermal_cov{i}=1/(J-1)*[En.thermal{i}-En.thermal_means{i}]*[En.thermal{i}-En.thermal_means{i}]';
    [U,D,V]=svd(En.thermal_cov{i});
    En.thermal_R{i}=U*sqrt(D);
end


%%% use posterior to make predictions on the interval [T,2T]

pred.int=Meas.HF.int(ind_m+1:2*Wall.Nt)';
pred.ext=Meas.HF.ext(ind_m+1:2*Wall.Nt)';
inv_sqrt_Gamma_pred.int=diag(1./(Meas.stda.int(ind_m+1:2*Wall.Nt)));
inv_sqrt_Gamma_pred.ext=diag(1./(Meas.stda.ext(ind_m+1:2*Wall.Nt)));
J_res=1e4; %Number of resample particles

for en=1:J_res
    %%Resample from Gaussian approximation
         for i=1:5
             thermal{i}=En.thermal_means{i}(:,1)+ En.thermal_R{i}*randn(length(En.thermal{i}(:,1)),1);
         end    
    HF=Heat_FEM(Grid.hx,Grid.Nx,2*Wall.Nt,Wall,exp(thermal{2}),...'
        exp(thermal{1}),thermal{3},exp(thermal{4}),exp(thermal{5}),Meas.T_int,Meas.T_ext,1);
    Heat_pred.int(:,en)=HF.int(ind_m+1:2*Wall.Nt);
    Heat_pred.ext(:,en)=HF.ext(ind_m+1:2*Wall.Nt);
    
end

%%Compute chi-square (see equation (27)) in De Simon et-al
Pos.Err_pred.int=norm(inv_sqrt_Gamma_pred.int*((mean(Heat_pred.int') - pred.int)'))^2;
Pos.Err_pred.ext=norm(inv_sqrt_Gamma_pred.ext*((mean(Heat_pred.ext') - pred.ext)'))^2;
l_int=prctile(Heat_pred.int',2.5);
l_ext=prctile(Heat_pred.ext',2.5);
u_int=prctile(Heat_pred.int',97.5);
u_ext=prctile(Heat_pred.ext',97.5);
%%Compute Average Interval Score (see equation (28)) in De Simon et-al;
aa=0.05;
Pos.Score.int=mean(u_int-l_int+2/aa*(l_int-pred.int).*(pred.int<l_int)+2/aa*(pred.int-u_int).*(u_int<pred.int));
Pos.Score.ext=mean(u_ext-l_ext+2/aa*(l_ext-pred.ext).*(pred.ext<l_ext)+2/aa*(pred.ext-u_ext).*(u_ext<pred.ext));
%%Compute posterior predictive credible intervals and predictive mean of
%%heat fluxes
Pos.Predictive.low=[prctile(Heat_pred.int',2.5);prctile(Heat_pred.ext',2.5)];
Pos.Predictive.high=[prctile(Heat_pred.int',97.5);prctile(Heat_pred.ext',97.5)];
Pos.Predictive.mean=[mean(Heat_pred.int');mean(Heat_pred.ext')];





function phi_r=find_temp(J_thresh,log_W,phi_r_1)
%%% function to compute tempering parameter defined by equation B.5 (using
%%% a bisection algorithm. 
%Inputs: 
%%% J_thresh: Effective sample size (see equation B.5)
%%% log_W: logarithm of weights (negative log-likelihood)
%%% phi_r_1: previous tempering parameter

%Outputs:
%%% phi_r: tempering parameter that satisfies B.5
phi_r=1;
W_new=exp(-log_W*(phi_r-phi_r_1));
%%% If weights are very small:
while (norm(W_new)<1e-10)
    phi_r=(phi_r+phi_r_1)/2;
    W_new=exp(-log_W*(phi_r-phi_r_1));
end
W_new=W_new/sum(W_new);
ESS=1/sum(W_new.^2);
phi_t=phi_r_1;
max_iter=1000;
iter=0;
%%%Bisection algorithm
if (ESS<J_thresh)
    while (abs(ESS-J_thresh)>1e-1)&&(iter<max_iter)
        iter=iter+1;
        if (iter==max_iter)
            disp('maximum number of iterations reached')
        end
        if (ESS<J_thresh)
            phi_r_old=phi_r;
            phi_r=(phi_t+phi_r)/2;
            W_new=exp(-log_W*(phi_r-phi_r_1));
            W_new=W_new/sum(W_new);
            ESS=1/sum(W_new.^2);
        else            
            phi_t=phi_r;
            phi_r=(phi_r_old+phi_r)/2;
            W_new=exp(-log_W*(phi_r-phi_r_1));
            W_new=W_new/sum(W_new);
            ESS=1/sum(W_new.^2);
        end
    end
end

