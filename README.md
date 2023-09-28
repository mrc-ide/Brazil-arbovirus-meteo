# Arbovirus_Brazil_reg

This repository contains code to recreate the spatiotemporal regression modelling study in the preprint 'Spatiotemporal relationships between climate and extreme weather events and arbovirus human infections across Brazil' (https://doi.org/________). 

## Repository overview
This repository contains scripts to analyse the following dataset (_____________________________). We have made the dataset publicly available under a ____ liscence, please cite both the paper and the dataset when including them in your work. The data contains the 
- run a univariate analysis
- build mulivariate models in a stepwise bottom-up approach
- fit multivariate models using the R-INLA package and assess how well

## Code
All analyses can be run from the main.R file, which calls scripts to run the analysis (contained in the scripts folder) and functions (contained in the R folder). The purpose of each R file is described below.

### Scripts
#### _simulate_data.R_
  
Simulate 540 serological datasets - each dataset is a linelist of individuals with age and simulated antibody titre level.
  
#### _simulated_data_thresholds.R_
  
Estimates optimal antibody titre thresholds per dataset to classify the titres as seropositive (indicating a previous infection) or seronegative (indictating no previous infection) which is necessary for the catalytic models.
  
#### _to_run_catalytic_models.R_

Main script to fit the catalytic models.

#### _plotting_catalytic_models_1.R_

Code to generate and save .csv files of the seroprevalence and FOI estimates, and to generate plots of the estimated age-specific seroprevalence.

#### _to_run_mixture_model.R_

Main script to fit the mixture model.

#### _plotting_mixture_model_1.R_

Code to generate plots of the mixture model fitted to the histogram of the simulated antibody titre distributions.

#### _plotting_mixture_model_2.R_

Code to generate and save .csv files of the seroprevalence and FOI estimates, and to generate plots of the estimated age-specific seroprevalence.

#### _explore_results_1.R_

Mixture model results. Assess the ability of the mixture model to correctly assign the distribution family to the seronegative and seropositive compartments of the simulated datasets. Compare and plot the estimated mean titre scores, distribution standard deviations and the FOI to the 'true' simulated values. Further, calculate uncertainty and bias in the FOI and seroprevalence estimates.  

#### _explore_results_2.R_

Catalytic model results. Calculate and plot uncertainty and bias in the FOI and seroprevalence estimates.

#### _explore_results_3.R_

Mixture and Catalytic model output comparisons. Compare the estimated FOI and seroprevalence values, bias and uncertainty between the three models and generate comparison plots.


### R
#### _simulating_data_functions.R_

Code for the simulate() function to simulate the data and for the est_serostatus() function which is minimised to calculate classification thresholds in the _simulated_data_thresholds.R_ script

#### _cat_model_functions.R_

Functions to process the categorised data and fit the time-varying FOI and time-constant FOI catalytic models.

#### _mix_fitting_functions.R_
 
Functions for fitting the mixture model to the titre distributions (step 1 of mixture model).

#### _mix_calc_functions.R_
