function [Rel, Pot, CRPS_tot] = CRPS_decomposition_hersbach(ens, obs)

%-------------------------------------------------------------------------------------------------------------
% Rachel Bazile, March 20 2018, from the original code of Marie-Amelie Boucher - September 2008 - Stage IREQ
%-------------------------------------------------------------------------------------------------------------
% This function decomposes the CRPS into reliability (Rel) and "potential" (Pot)
% components according to Hersbach (2000). The potential CRPS
% represents the best possible CRPS value that could be achieved, if the forecasts
% were perfectly reliable. The total (CRPS_tot) of Rel+Pot is the empirical CRPS according to the definition of
% Hersbach (2000)
%
% ------------------------------------------------------------------------------------------------------------ 
%   INPUT
%       ens :  matrix of size (m x n) with m the number of validity dates and n the number of members
%       obs  : vector of size (m x 1) containing the observations matching the forecasts
%
%   OUTPUT
%       Rel :  Scalar - the reliability component of the CRPS
%       Pot  : Scalar - the potential component of the CRPS
%       CRPS_tot :  Scalar - total (empirical) CRPS following the definition of Hersbach (2000). 
%
%
% Hersbach, H., 2000. Decomposition of the continuous ranked probability score for ensemble 
% prediction systems. Weather Forecast. 15, 550Ð570.
%------------------------------------------------------------------------------------------------------------

[M,N] = size(ens);
alpha = zeros(M,N);
beta = zeros(M,N);

for i = 1 : M  % The loop is on the number of forecast-observation pairs in the available archive 
    
    if ~isnan(obs(i)) % If the observation in NaN, it is not possible to compute alpha and beta 
        
        ens_sort = sort(ens(i,:));
                 
        for k = 1 : N+1   % The loop is on the number of members + 1 (the number of intervals in Hersbach's empirical distribution)
                      
            if k == 1
                if obs(i) < ens_sort(1)
                    alpha(i,k) = 0;
                    beta(i,k) = ens_sort(1) - obs(i);
                else
                    alpha(i,k) = 0;
                    beta(i,k) = 0;
                end
            elseif k == N+1
                if obs(i) > ens_sort(N)
                    alpha(i,k) = obs(i) - ens_sort(N);
                    beta(i,k) = 0;
                else
                    alpha(i,k) = 0;
                    beta(i,k) = 0;
                end 
            else
                if obs(i) > ens_sort(k)
                    alpha(i,k) = ens_sort(k) - ens_sort(k-1);
                    beta(i,k) = 0;
                elseif obs(i) < ens_sort(k-1)
                    alpha(i,k) = 0;
                    beta(i,k) = ens_sort(k) - ens_sort(k-1);
                elseif ((obs(i) >= ens_sort(k-1)) && (obs(i) <= ens_sort(k)))
                    alpha(i,k) = obs(i) - ens_sort(k-1);
                    beta(i,k) = ens_sort(k) - obs(i);
                else
                    alpha(i,k) = NaN;
                    beta(i,k) = NaN;
                end  
            end    
        end
        
    else
        alpha(i,:) = NaN;
        beta(i,:) = NaN;
    end
end

alpha1 =nanmean(alpha,1);
beta1 = nanmean(beta,1);

g = alpha1+beta1;  
o = beta1./(alpha1+beta1);

Rel = nansum(g.*(o-(0:N)/N).^2);
Pot = nansum(g.*o.*(1-o));
CRPS_tot = Rel+Pot;

%
