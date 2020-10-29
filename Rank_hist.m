function [Rel_freq, bins]=Rank_hist(ens, obs)

%Marie-Amelie Boucher
%Mars 2006
%Universite Laval
%----------------------------------------------------------------------------------------------------------
% This function computes the required variables to plot rank histograms (Hamill and Colucci, 1997; Talagrand et al., 1997)
%
% ********. NOTE: The function does *not* plot the histogram! Should you want to plot it, 
% you would have to type the command "bar(bins, Rel_freq)", or just "bar(Rel_freq)".
%
% Hamill, T.M. and S.J. Colucci. 1997. "Verification of Eta-RSM short-range ensemble forecasts." Monthly Weather Review, 125, 1312-1327
% Talagrand, O., R. Vautard and B. Strauss. 1997. "Evaluation of probabilistic prediction systems." Proceedings of the ECMWF Workshop on
% predictability, ECMWF, Reading, UK, 1-25
%
% INPUTS: 
%
%--------------> ens    =   a (n x m) matrix, where n is the number of time steps (validity dates) and m is the ensemble size
%
%--------------> obs    =   a (n x 1) vector, with n the number of time steps (validity dates).
%
% OUTPUTS:
% 
%--------------> Rel_freq     =   Relative frequency for each bin (frequency with which the observation falls into each bin)
%--------------> bins         =   One bin per ensemble member, plus one.

%-----------------------------------------------------------------------------------------------------------

ens=double(ens);

% Remove NaNs
indnan=find(~isnan(ens(:,1)));

ens_NoNaN=ens(indnan,:);
obs_NoNaN=obs(indnan,:);

[m,n]=size(ens_NoNaN);

percentiles = NaN(m,1);


for i=1:m
    
    augmented_ens = [ens_NoNaN(i,:),obs_NoNaN(i)];
    val = obs_NoNaN(i);

    if ~isnan(val)   % If the observation is NaN, it is not possible to rank it.

        sorted_augmented_ens = sort(augmented_ens);

        index = find(sorted_augmented_ens==val) ;    %This provides the rank occupied by the observation
        [~,b] = size(index);

        if b>1  % In case of rank equality (i.e., two or more values of the augmented ensemble are perfectly equal), we perform a random draw to determine the rank of the observation
            tirage = ceil((rand)*b);
            index = index(tirage);
        end

        ranks(i) = index;
    else
        ranks(i) = NaN;
    end

end


% Computing the relative frequency for each rank.

[N,x] = hist(ranks,(n+1));
bins = 1:n+1;

Rel_freq = N/m;
