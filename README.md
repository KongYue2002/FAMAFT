We have attached one R file and two folders to reproduce all the results in the simulations and real data analysis of the manuscript "Factor-augmented Multivariate Accelerated Failure Time Model for Multiple Events Data". The details are listed as follows:


File:
    
functions.R 
— This script contains all major functions used in the numerical analysis of the manuscript, including functions for random data generation and algorithm implementation. This file should be sourced when running the main R scripts


Folders:

(1) simulation 
— This folder contains five files used to reproduce all simulation results, including tables and figures. 

sim_code_main.R:	The main simulation script. It calls functions.R and contains the core loop for data generation and parameter estimation.

beta_hat_results.R: The R script for organizing and saving the simulation results for estimated regression coefficients.

factor_and_loading_results.R: The R script for saving the simulation results for estimating factors and loadings, and selecting factors number based on eigenvalue ratio.

prediction_results.R: The R script for saving the simulation results for predicting the true survival times.

time_results.R: The R script for saving the computational times (in seconds) for different methods.

(2) real_data_analysis
— This folder contains one folder and two R files used to reproduce all the results in the real data  analysis of the of the manuscript.

data: The folder that contains the input data files used for the real-data analysis.

real_code_main.R:	The main script for the real-data analysis. It loads the data and runs the estimation procedures.

results.R:	A separate script for organizing, summarizing, and saving the output results (including the estimated regression coefficients, factors and loadings, C-indices and computation times).

