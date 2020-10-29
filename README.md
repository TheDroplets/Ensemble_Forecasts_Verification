# Ensemble_Forecasts_Verification
Verification scores and tools for ensemble forecasts
Created on Thu Oct 22 12:47:06 2020

Marie-Amelie Boucher, marie-amelie.boucher@usherbrooke.ca

This library contains the following functions:

- CRPS

- CRPS_decomposition_Hersbach

- Log_score

- Rank_hist

- Reliability

They were originally coded in Matlab by me, mainly during my MSc and PhD studies. Throughout the years, students have greatly helped to improve those functions (mostly but not exclusively Rachel Bazile). 

IMPORTANT NOTE: NONE of those functions actually *plot* something. It does, however, computes elements that can be plotted (e.g. the variables needed to plot a rank histogram)


## Comments on specific functions:

### CRPS

This function computes the CRPS for a series of ensemble forecasts (ens) and corresponding observations (obs). The user can choose which distribution to fit to the ensemble for 
computation.

INPUTS

- ens: a (n x m) matrix, where n is the number of time steps (validity dates) and m is the ensemble size

- obs: a (n x 1) vector, with n the number of time steps (validity dates). 

- distribution = 'emp': an empirical distribution is computed from the ensemble forecasts

- distribution = 'normal_exact':  a normal (Gaussian) distribution is fitted to the ensemble forecasts. The CRPS is computed using an exact formulation derived from the equation of the normal distribution. 

or 

distribution = 'gamma_approx': a gamma distribution is fitted to the ensemble forecasts. The CRPS is computed using an approximate formulation as per Gneiting and Raftery (2007) equation (21)

Tilmann Gneiting and Adrian E. Raftery (2007): Strictly Proper Scoring Rules, Prediction, and Estimation, Journal of the American Statistical Association, 102:477, 359-378

OUTPUT:

- CRPS: a scalar, representing the mean CRPS for the time series of forecasts 

Original by Marie-Amelie Boucher, July 2005, Universite Laval
Modified by Rachel Bazile and Luc Perreault in June 2017: 
The computation for the empirical version of the CRPS cannot be based on just a few
ensemble members so we need to "subsample" to have enough members



## CRPS decomposition hersbach

This function decomposes the CRPS into reliability and "potential" 
components according to Hersbach (2000). The potential CRPS represents the best possible CRPS value that could be achieved, if the forecasts were perfectly reliable. The total (CRPS_tot) is the empirical CRPS according to the definition of Hersbach (2000)

Hersbach, H., 2000. Decomposition of the continuous ranked probability score for ensemble 
prediction systems. Weather Forecast. 15, 550–570.

inputs: 

- ens:    mxn matrix; m = number of records (validity dates)  n = number of members in ensemble 

- obs:    mx1 vector; m = number of records (validity dates, matching the ens) 

outputs:   

- CRPS_tot:       Scalar -> reliability + potential 

- reliability:    Scalar -> the reliability component of the CRPS

- potential:      Scalar -> the potential CRPS that would be reached for a perfectly reliable system


## Log_score

This function computes the logarithmic (or ignorance) score. Predictive distributions can be considered as Gaussian, Gamma distributed, Empirical or "Loi des fuites"
(a Gamma distribution + a Dirac at zero, suitable for daily precip), and Kernel distribution.

inputs: 

- ens:    mxn matrix; m = number of records (validity dates)  n = number of members in ensemble 

- obs:    mx1 vector; m = number of records (validity dates, matching the ens)

case:           

- 'Normal'

- 'Gamma'

- 'Kernel'

- 'Fuites'  is made for daily precipitation exclusively 

- 'Emp'

- thres:          probability density threshold below which we consider that the event was missed by the forecasting system. This value must be small (e.g.: 0.0001 means that f(obs) given the forecasts is only 0.0001 --> not forecasted). 

opt_case         

- if 'case' = 'Fuites', opt_cas is the threshold to determine data which contributed to gamma distribution and those who are part of the Dirac impulsion

- if 'case' = 'Emp', opt_cas needed is the number of bins in which to divide the ensemble, by default, it will be the number of members (Nan excluded). opt_cas have to be an integer superior to 1.

outputs:

- S_LOG:           the logarithmic score (scalar)

- ind_miss:        Boleans to point out days for which the event was missed according


Reference:

'Emp' case is based on Roulston and Smith (2002) with modifications -> quantile and members with similar values

History

MAB June 19: Added 2 cases for the empirical distribution: the observation can either be the smallest or the largest member of the augmented ensemble, in which case we can't use the "DeltaX = X(S+1)-X(S-1);" equation.


## Rank_hist

This function computes the required variables to plot rank histograms 
(Hamill and Colucci, 1997; Talagrand et al., 1997) 
********. NOTE: The function does *not* plot the histogram! It just computes what you need to plot it

Hamill, T.M. and S.J. Colucci. 1997. “Verification of Eta-RSM short-range ensemble forecasts.” Monthly Weather Review, 125, 1312-1327

Talagrand, O., R. Vautard and B. Strauss. 1997. “Evaluation of probabilistic prediction systems.” Proceedings of the ECMWF Workshop on predictability, ECMWF, Reading, UK, 1-25

inputs: 

- ens    =   a (n x m) matrix, where n is the number of time steps (validity dates) and m is the ensemble size

- obs    =   a (n x 1) vector, with n the number of time steps (validity dates).

outputs:  
 
- Freq     =   Relative frequency for each bin (frequency with  which the observation falls into each bin)

- bins     =   One bin per ensemble member, plus one.


## Reliability

This function computes what is needed to plot (modified) reliability diagrams.
Note that this function does *not* adopt the definition of Murphy and Winkler (1977), 
which would require the transformation of the ensemble into a binary forecast using 
a threshold.

Instead, we verify whether each confidence interval (i.e. 10%,..., 90%) corresponds
to its definition: for instance the 10% confidence interval should include 10% of 
observations.

This function *does not* plot the diagram. It just computes what you need to plot it.

Finally, as an additional information, we also compute the mean length of each confidence interval. A variant of the reliability diagram could then be: "plot(length, effective_coverage)", as in Boucher et al. (2010)

References:

Murphy AH and Winkler RL (1977) Reliability of Subjective Probability Forecasts of Precipitation and Temperature, Journal of the Royal Statistical Society: Series C (Applied Statistics) 26: 41-47

Boucher M-A, Laliberte J-P and Anctil F. (2010) An experiment on the evolution of an ensemble of neural networks for streamflow forecasting, Hydrol. Earth Syst. Sci., 14, 603–612

inputs: 

- ens:    mxn matrix; m = number of records (validity dates)  n = number of members in ensemble 

- obs:    mx1 vector; m = number of records (validity dates, matching the ens) 

outputs:

- nominal_coverage    : nominal probability of the intervals

- effective_coverage  : effective coverage of the intervals

- Length              : mean effective length of the intervals


