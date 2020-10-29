function [S_LOG, ind_miss]=Log_score(ens, obs, distribution, thres, options)

% Coded by Marie-Amelie Boucher and Rachel Bazile (February 2017)
%--------------------------------------------------------------------------
% This function computes the logarithmic (or ignorance) score. Predictive distributions can
% be considered as Gaussian, Gamma distributed, Empirical, Kernel distribution (also empirical) or "Loi des fuites"
% (a Gamma distribution + a Dirac at zero, suitable for daily precipitation forecasts).
%
% INPUTS:
%       > ens : n * m matrix of ensemble forecasts. n is the number of
%         days (validity dates) and m the number of members.
%       > obs : n * 1 vector of observations matching the forecasts
%       > distribution : 'Normal' or 'Gamma' or 'Kernel' or 'Fuites' or 'Emp'
%              'Fuites' is mostly suitable for daily precipitation forecasts
%       > thres : probability density threshold below which we consider that the
%         event was missed by the forecasting system. This value must be
%         small. For the 'Emp' (empirical) distribution, the threshold is replaced by minimal and
%         maximal bounds (vector of 2x1 elements).
%         By default, thres = 0.00001 and the logarithmic score is unbounded. 
%       > options : options relative to specific distributions:
%                   >   if 'distribution' = 'Kernel', 'options' is a very small
%                       value, needed to comprise the ensemble between
%                       min(ens_prev)-value and max(ens_prev) + value
%                   >   if 'distribution' = 'Fuites', 'options' is a
%                       threshold (precipitation) value below which the ensemble members
%			are treated as zeros. The purpose is to define which portion of 
%			the ensemble is used to fit the gamma distribution
%                   >   if 'distribution' = 'Empirical', 'options' is the
%                       number of bins in which to divide the ensemble. By
%                       default, it will be the number of members (NaN
%                       excluded). 'options' has to be an integer larger than 1
%                       
%
% OUTPUTS:
%       > S_LOG: the logarithmic score (scalar)
%       > ind_miss: the number of missed events, i.e. the number of times the observation 
%		falls completely outside of the ensemble (scalar)
%
% Reference:
%   The 'Emp' case is based on Roulston and Smith (2002) with
%   modifications to discriminate between quantile and members with similar values
%
% MARK S. ROULSTON AND LEONARD A. SMITH (2002) "Evaluating Probabilistic Forecasts 
% Using Information Theory", Monthly Weather Review, 130, 1653-1660.
%
%%--------------------------------------------------------------------------
% History
%
% MAB June 2019: Added 2 cases for the empirical distribution: the
% observation can either be the smallest or the largest member of the
% augmented ensemble, in which case we can't use the "DeltaX = X(S+1) -
% X(S-1);" equation.
% --------------------------------------------------------------------------
narginchk(3, 5)

[n,~]=size(ens);
loga=nan(n,1);
ind_miss=nan(n,1);

if length(obs) ~= n
    error('The length of the record of observations doesn''t match the length of the forecasting period')
end

if ~strcmp(distribution,'Emp') && (~exist('thres','var') || thres == 0 )
    thres = 0.00001;
    disp('The logarithmic score is unbounded')
elseif ~strcmp(distribution,'Emp') && ((thres <= 0) || (thres > 1))
    error('The threshold has to be positive and smaller than 1.')
end

switch distribution
    
    case {'Emp'}
        
        if ~exist('options','var') % If a number of bins is specified
            disp('Bins used for the empirical method are those defined by the ensemble members')
        elseif options < 2 || (round(options) ~= options) % This doesn't work well
            error('Invalid format for options.')
        end
        
        if ~exist('thres','var')  || (length(thres) ~= 1)
            error('Invalid format for the thresholds. For the Emp case, thresholds is a 2x1 vector.')
        end
        
        for j = 1 : n
            
            if sum(isnan(ens(j,:))) < length(ens(j,:))
                
                if (min(ens(j,:)) <= obs(j)) && (obs(j) <= max(ens(j,:)))
                    ind_miss(j) = 0;
                else
                    ind_miss(j) = 1;
                end
                
                % Suppress NaN from the ensemble to determine the number of
                % members
                prev = sort(double(ens(j,~isnan(ens(j,:)))));
                
                if exist('options','var') % If a number of bins is specified
                    prev = quantile(prev,0:1/options:1);
                end
                
                N = length(prev) ;
                
                if length(unique([prev,obs(j)])) == 1 % if all the members and the observation are identical -> perfect forecast
                    proba_obs = 1;
                    S = [];
                else
                    if length(unique(prev)) ~= length(prev) % if some members are identical, modify their values very slightly 
                        rng('shuffle','v5normal' )
                        uni_prev = unique(prev);
                        [ind1,ind2] = histc(sort(prev),unique(sort(prev)));
                        ind = find(ind1 > 1);
                        for k = ind
                            new_val = uni_prev(k) + 0.01.*rand(ind1(k)-1,1);
                            prev(ind2 == k) = [uni_prev(k) ; new_val];
                        end
                    end
                    X = sort([thres,prev,obs(j)]);
                    S = find(X == obs(j));
                end
                
                if length(S) == 1 % if the observation falls between two members OR occupies the first or last rank
                    
                    if S==length(X) % If the observation is the largest member of the augmented ensemble
                    proba_obs = thres;   
                        
                    elseif S==1  % If the observation is the smallest member of the augmented ensemble
                    proba_obs = thres;  
                    
                    else
                        
                        DeltaX = X(S+1) - X(S-1);
                        proba_obs = min(1/(DeltaX*(N+1)),1);
                    end
                    
                elseif length(S) == 2  % if the observation is identical to one member, choose the maximum probability density between the two (observation or member)
                    DeltaX1 = X(S(2)+1) - X(S(2));
                    DeltaX2 = X(S(1)) - X(S(1)-1);
                    DeltaX = min(DeltaX1,DeltaX2);
                    proba_obs = min(1/(DeltaX*(N+1)),1);
                end
                loga(j) = -log2(proba_obs);
                
            else
                
                loga(j) = NaN ;
                ind_miss(j) = NaN ;
            end
        end
        
        
        
    case {'Normal'}
        
        for j = 1 : n
            
            if sum(isnan(ens(j,:))) < length(ens(j,:))
                pd = fitdist(ens(j,:)',distribution);
                proba_obs = min(pdf(pd,obs(j)),1);
                
                if proba_obs >= thres
                    ind_miss(j) = 0;
                    loga(j)=-log2(proba_obs);
                else
                    loga(j) = -log2(thres);
                    ind_miss(j) = 1;
                end
                
            else
                
                loga(j) = NaN ;
                ind_miss(j) = NaN;
                
            end
            
        end
        
    case {'Gamma'}
        
        if sum(sum(ens <= 0)) == 0
            
            for j = 1 : n
                
                if sum(isnan(ens(j,:))) < length(ens(j,:))

                    pd = fitdist(ens(j,:)',distribution);
                    proba_obs = min(pdf(pd,obs(j)),1);
                    
                    if proba_obs >= thres
                        ind_miss(j) = 0;
                        loga(j)=-log2(proba_obs);
                    else
                        loga(j) = -log2(thres);
                        ind_miss(j) = 1;
                    end
                    
                else
                    
                    loga(j) = NaN ;
                    ind_miss(j) = NaN;
                    
                end
                
            end
            
        else
            error('Forecasts contain zeros. You must choose a different distribution.')
        end
        
    case {'Kernel'}
        
        if ~exist('options','var') || (options == 0)
            support = 'unbounded';
            disp('The Kernel distribution is unbounded.')
        end
        
        for j = 1 : n
            
            if sum(isnan(ens(j,:))) < length(ens(j,:))
                
                support = [min(ens(j,:)) - options , max(ens(j,:)) + options];
                pd = fitdist(ens(j,:)',distribution,'support',support);
                proba_obs = min(pdf(pd,obs(j)),1);
                
                if proba_obs >= thres
                    ind_miss(j) = 0;
                    loga(j)=-log2(proba_obs);
                else
                    loga(j) = -log2(thres);
                    ind_miss(j) = 1;
                end
                
            else
                
                loga(j) = NaN ;
                ind_miss(j) = NaN ;
            end
            
        end
        
    case 'Fuites'
        
        if ~exist('options','var')
            error('Option missing for ''Fuites'' distribution.')
        end
        
        for j=1:n
            if sum(isnan(ens(j,:))) < length(ens(j,:))
                ind_non_null = find(ens(j,:) > options);
                prop_null = (length(ens(j,:))-length(ind_non_null)) / length(ens(j,:));
                if obs(j) <= options
                    proba_obs = prop_null;
                else
                    ens_non_null = ens(j,ind_non_null);
                    ab = gamfit(double(ens_non_null));    % Fitting gamma parameters (max. likelihood method))
                    proba_obs = min(gampdf(obs(j),ab(1),ab(2))*(1-prop_null),1);
                end
                if proba_obs >= thres
                    ind_miss(j) = 0;
                    loga(j)=-log2(proba_obs);
                else
                    loga(j) = -log2(thres);
                    ind_miss(j) = 1;
                end
                
            else
                
                loga(j) = NaN ;
                ind_miss(j) = NaN ;
            end
            
        end
        
    otherwise
        
        error('Choice of distribution type in ''distribution'' is incorrect. Possible options are : ''Normal'',''Gamma'',''Kernel'', ''Emp'' or ''Fuites''')
        
end

S_LOG = nanmean(loga);
ind_miss = nansum(ind_miss) ; 

%
