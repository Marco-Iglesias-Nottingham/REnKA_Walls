# REnKA_Walls
MATLAB code for the sequential estimation of the thermophysical properties of a wall given measurements of surface heat flux and near air temperatures. The code uses the ensemble Kalman inversion algorithm, as outlined in the paper "Quantifying uncertainty in thermophysical properties of walls by means of Bayesian inversion" by De Simon, Iglesias, Jones, and Wood (To appear in Energy and Buildings). The code is provided in the context of the synthetic experiment described in Sections 3.1 and 4.1-4.2.

The following files are provided:

(1) "generate_synthetic_data"
	Generate synthetic data as described in Section 3.1.
	
(2) "generate_prior_ensemble"
	Generates an ensemble of the prior as discussed in Appendix A.
	
(3) "REnKA"
	Takes the prior ensemble and synthetic measurements to compute approximate approximate posteriors of thermophysical properties as well as predictive distributions of heat flux. This function is based on the ensemble Kalman inversion algorithm described in Appendix B
		
(4)  "visualise_results"
	Visualises synthetic data from (1) as well as some statistics from the posterior computed via (3). These include (i) credible intervals and mean of (k,c,T_{0}) computed at a given assimilation time,  (ii) running mean and credible intervals of R_I,R_E, U-value and C-value, (iii) prior and posterior credible interval of the predictive distributions of surface heat flux.

(5) "cova_matrix"
	Code to compute the Karhunen-Loeve decomposition of the prior covariance matrices of k, c and T_{0}. This code is used in (2) to generate the prior ensemble of (k,c,T_{0}) as well as in (1) to produce a stochastic perturbation to the synthetic near air temperatures. 

(6) "Heat_FEM"
	Finite Element solver of the 1D heat equation used to generate synthetic data. This solver is also used to compute the forward model evaluations within the ensemble Kalman inversion framework.
	
(7) "Driver.m"
	Main script that demonstrates the use of (1)-(7) and to reproduce some of the results from sections 4.1-4.2.

