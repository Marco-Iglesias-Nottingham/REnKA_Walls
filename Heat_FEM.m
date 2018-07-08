function out=Heat_FEM(hx,Nx,Nt,Wall,K,c,P0,R_int,R_ext,T_int,T_ext,option)
%%% Finite element solver for the 1D heat equation (see details in section 2.5
%%% of "Quantifying uncertainty in thermophysical properties of walls 
%%% by means of Bayesian inversion" by De Simon, Iglesias, Jones, 
%%% and Wood (To appear in Energy and Buildings)

%Inputs: 
%%% hx: element size (assumes regular  
%%% Nx: Number of elements. 
%%% Nt: final time step. 
%%% Wall: structure with Wall's geometry and time-step (Wall.tau)
%%% K: thermal conductivity
%%% c: volumetric heat capacitance
%%% P0: Initial temperature distribution.
%%% R_int: internal surface resistance
%%% R_ext: external surface resistance
%%% T_int: internal temperature (at each time step)
%%% T_ext: external temperature (at each time step)
%%% option: 0 or 1

%Outputs:
%%% out: If option 0 then out is the final temperature distribution after Nt time-steps.
%%% If option 1 then out is the internal and external surface heat fluxes at all 
%%% time-steps between 0 and Nt.

%%%Note: This FEM solver uses linear basis functions for the temperature
%%%profile and piecewise constant for k and c. These assumptions have been
%%%hard-coded for the sake of efficiency in the generation of the stiffness
%%%and mass matrices. A more general solver can be used to replace Heat_FEM.




rhs=zeros(Nx+1,1);
K_hx(1)=1/hx/R_int;
K_hx(2:Nx+1)=K/hx^2;
K_hx(Nx+2)=1/hx/R_ext;
x1=reshape(K_hx(1:Nx+1),Nx+1,1); x2=reshape(K_hx(2:Nx+2),Nx+1,1);
DiagVecs=[-x2,x1+x2,-x1];
DiagIndx=[-1,0,1];
A=spdiags(DiagVecs,DiagIndx,Nx+1,Nx+1);


c_h(1)=0;
c_h(2:Nx+1)=c;
c_h(Nx+2)=0;
c1=reshape(c_h(1:Nx+1),Nx+1,1); c2=reshape(c_h(2:Nx+2),Nx+1,1);
DiagVecs2=[c1/6,(c1+c2)/3,c2/6];
DiagIndx2=[-1,0,1];
C=spdiags(DiagVecs2,DiagIndx2,Nx+1,Nx+1);


P=P0;
V_int(1)=(T_int(1)-P(1))*K_hx(1)*hx;
V_ext(1)=(P(Nx+1)-T_ext(1))*K_hx(Nx+2)*hx;
v1=K_hx(1)*T_int;
v2=K_hx(Nx+2)*T_ext;
for t=1:Nt-1 
    rhs(1)=v1(t+1);
    rhs(Nx+1)=v2(t+1);
    rhs2=rhs+C/Wall.tau*P;
    P=(A+C/Wall.tau)\rhs2;  
    V_int(t+1)=(T_int(t+1)-P(1))*K_hx(1)*hx;
    V_ext(t+1)=(P(Nx+1)-T_ext(t+1))*K_hx(Nx+2)*hx;
end

HF.int=V_int';
HF.ext=V_ext';

if (option==1)
    out=HF;    
elseif (option==0)
    out=P;
end


