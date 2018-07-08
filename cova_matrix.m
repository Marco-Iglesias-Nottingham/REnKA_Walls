function R=cova_matrix(x,prior)
%%% Generate the Karhunen-Loeve (KL) decomposition of the prior covariance 
%%% with Whittle-Matern correlation function with parameters sigma, nu and L
 
%Inputs: 
%%%% x: locations of centres of cells/nodes
%%% prior: structure with hyperparameters sigma, nu and L of a 
%%% Whittle-Matern prior covariance matrix C 
%Outputs:
%%% R: KL decomposition matrix of C
%%% Use xi=R*randn(length(x),1) to generate a
%%% Gaussian random field with zero mean and covariance C.

sigma=prior(1);
nu=prior(2);
L=prior(3);
N=length(x);

C=zeros(N,N);
%%% construct the prior error covariance
for i=1:N
    h=abs((x(i)-x(1:N)));
    C(i,1:N)=  sigma^2*2^(1-nu)/gamma(nu)*(h(1:N)/L).^nu.*besselk(nu,h(1:N)/L);
end
for i=1:N
    C(i,i)=  sigma^2;
end


[U,D,V]=svd(C);
%%Compute KL decomposition matrix
R=U*sqrt(D);

