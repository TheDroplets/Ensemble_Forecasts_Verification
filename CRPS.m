function [S_CRPS]=CRPS(ens, obs, distribution)

%--------------------------------------------------------------------------
% This function computes the CRPS for a series of ensemble forecasts (ens) and corresponding
% observations (obs). The user can choose which distribution to fit to the ensemble for 
% computation.
%
% INPUTS
%-----------------------------------------------------------------------------------
% ens: a (n x m) matrix, where n is the number of time steps (validity dates) and m is the ensemble size
% obs: a (n x 1) vector, with n the number of time steps (validity dates). 
% distribution = 'emp': an empirical distribution is computed from the ensemble forecasts
% distribution = 'normal_exact':  a normal (Gaussian) distribution is fitted to the ensemble forecasts. 
% The CRPS is computed using an exact formulation derived from the equation of the normal distribution. 
% or 
% distribution = 'gamma_approx': a gamma distribution is fitted to the ensemble forecasts. The CRPS
% is computed using an approximate formulation as per Gneiting and Raftery (2007) equation (21)
% Tilmann Gneiting and Adrian E. Raftery (2007): Strictly Proper Scoring Rules, Prediction, 
% and Estimation, Journal of the American Statistical Association, 102:477, 359-378
%
%-----------------------------------------------------------------------------------
% OUTPUT:
% S_CRPS: a scalar, representing the mean CRPS for the time series of forecasts 
%
% --------------------------------------------------------------------------------
% Original by Marie-Amelie Boucher, July 2005, Universite Laval
% Modified by Rachel Bazile and Luc Perreault in June 2017: 
% The computation for the empirical version of the CRPS cannot be based on just a few
% ensemble members so we need to "subsample" to have enough members
%--------------------------------------------------------------------------


% initialization

[m,~]=size(ens);
CRPS=nan(m,1);

switch distribution
    
    case 'emp' % Non-parametric estimation based on the empirical cumulative distribution of the ensemble. According to Luc Perreault's idea
        
        for j=1:m
            if ( sum(isnan(ens(j,:)))==0 && isnan(obs(j))==0 )
                ssample = sort(ens(j,:));
                stepsize = 1/length(ens(j,:));
                area1 = 0;
                area2 = 0;
                subsample1 = ssample(ssample <= obs(j));
                subsample2 = ssample(ssample > obs(j));
                n1 = length(subsample1);
                n2 = length(subsample2);
                for i=2:1:n1
                    area1 = area1 + ((i-1)*stepsize)^2*(subsample1(i)-subsample1(i-1));
                    if i==n1
                        area1 = area1 + (n1*stepsize)^2*(obs(j)-subsample1(i));
                    end
                end
                
                for i=2:1:n2
                    if i==2
                        area2 = area2 + (subsample2(i-1)-obs(j))*(stepsize*n2)^2;
                    end
                    area2 = area2 + ((n2-i+1)*stepsize)^2*(subsample2(i)-subsample2(i-1));
                    
                end
                CRPS(j) = (area1 + area2);
            else
                CRPS(j) = NaN;
            end
        end
        
        
    case 'normal_exact'
        [muhat,sigmahat] = normfit(ens');
        
        % Initialization
        
        vcr=nan(m,1);
        phi=nan(m,1);
        PHI=nan(m,1);
        
        for j=1:m
            vcr(j)=(obs(j)-muhat(j))/sigmahat(j);
            phi(j)=normpdf(vcr(j),0,1);
            PHI(j)=normcdf(vcr(j),0,1);
            
            CRPS(j)=abs(sigmahat(j)*  (   (1/sqrt(pi))-2*phi(j) - (vcr(j).*(2*PHI(j)-1))    ));
        end
        

    case 'gamma_approx'
        
        for j=1:m               %------------->It is not possible to fit a gamma distribution if there are negative values in the ensemble.
            a=find(ens(j,:)<=0);              
                                                                      
            if isempty(a)
                ens(j,:)=ens(j,:);
            else
                b=0.001;        % An arbitrary small value to replace negative values. This could be changed
                ens(j,a)=b;
            end
        end
        
        N=1000; % This is the number of draws from the fitted gamma distribution. This can be changed
       
        param_pdf_g_ajustee=zeros(m,2);
        
        for j=1:m
            phat= gamfit(ens(j,:)');
            param_pdf_g_ajustee(j,:)=phat;
            X=gamrnd(phat(1),phat(2),N,1);
            Xprime=gamrnd(phat(1),phat(2),N,1);
            CRPS(j)=nanmean(abs(X-obs(j)))-0.5*nanmean(abs(X-Xprime));
        end
               
end


S_CRPS=nanmean(CRPS);



