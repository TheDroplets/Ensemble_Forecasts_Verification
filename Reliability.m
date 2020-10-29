function [nominal_coverage, effective_coverage, length] = Reliability(ens, obs )

% Marie-Amelie Boucher, UQAC Oct. 2016
% -------------------------------------------------------------------------------
% This function computes what is needed to plot (modified) reliability diagrams.
% Note that this function does *not* adopt the definition of Murphy and Winkler (1977), 
% which would require the transformation of the ensemble into a binary forecast using 
% a threshold.
% Instead, we verify whether each confidence interval (i.e. 10%,..., 90%) corresponds
% to its definition: for instance the 10% confidence interval should include 10% of 
% observations.
%
% This function *does not* plot the diagram. To plot the diagram, one has to type 
% "plot(nominal_coverage, effective_coverage)"
%
% Finally, as an additional information, we also compute the mean length of each 
% confidence interval. A variant of the reliability diagram could be:
% "plot(length, effective_coverage)", as in Boucher et al. (2010)
%
% References:
%
% Murphy AH and Winkler RL (1977) Reliability of Subjective Probability Forecasts of 
% Precipitation and Temperature, Journal of the Royal Statistical Society: Series C
% (Applied Statistics) 26: 41-47
%
% Boucher M-A, Laliberte J-P and Anctil F. (2010) An experiment on the evolution of 
% an ensemble of neural networks for streamflow forecasting, 
% Hydrol. Earth Syst. Sci., 14, 603–612
%
%-------------------------------------------------------------------------------------
%
% INPUTS:
%   ens     : Ensemble forecasts (matrix of validity dates x nb of members)
%   obs     : Observations (vector of validity dates x 1)
%
%
% OUTPUTS:
%   nominal_coverage    : nominal probability of the intervals
%   effective_coverage  : effective coverage of the intervals
%   length              : mean effective length of the intervals
%-------------------------------------------------------------------------------------

% Nominal probability of the confidence intervals

p    = [0.05:0.05:0.45 0.55:0.05:0.95] ;     % probability bins
N   = length( p ) / 2 ;
nominal_coverage = diff( [p(1:N); p( end:-1:end-N+1 )] )' ;


[L,C] = size( ens ) ;
len   = zeros( L,N ) ;
eff   = zeros( L,N ) ;

for i = 1 : L
    
    q    = quantile( ens(i,:), p ) ;
    qmed = median( ens(i,:) ) ;
        
    % Compute lengths of intervals 
    
    len(i,:) = diff( [q(1:N); q( end:-1:end-N+1 )] ) ;
    
    % Locate the observation and associated percentile
    
    if obs(i) <= qmed
        eff(i,:) = q( 1:N ) <= obs(i) ;
    else
        eff(i,:) = q( end:-1:N+1 ) >= obs(i) ;
    end
end

% Compute averages

effective_coverage = nanmean( eff )' ;
length = nanmean( len )' ;



%
