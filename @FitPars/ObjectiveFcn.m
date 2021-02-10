function out = ObjectiveFcn(parms, this, ydata, sigma)
%% Calculate the cost (objective) function for data fitting
%

% -------------------------------------------------------------------------------------------------------------------------

    % check inputs
    if isempty(sigma) 
        if ~ismember(this.fitopts.noiseModel,{'none','default'}), error('%s: The input noise sigma cannot be empty',mfilename) ; end
    elseif ~isscalar(sigma) , error('%s: The input noise sigma should be a scalar here',mfilename) ; 
    end
    
    % calculate theory predicted signal
    E = feval(this.model.FcnSignal, parms, this.model) ; 
       
    % calculate cost function based on different noise models
    switch lower(this.fitopts.noiseModel)
        case {'none','default'} % no noise model
            out = ydata - E ; 
        case {'simple'} % simple rician noise model
            out = ydata - sqrt(E.^2 + sigma.^2) ; 
        case {'rician'} % Rician noise model
            out = - LogLikelihood(ydata, E, sigma) ; 
        case {'meanrician'} % adjust mean Rician noise
            out = ydata - sqrt(pi/2)*sigma.*Lhalf(-E.^2/2./sigma.^2) ;         
    end
    
    % out as a vector (lsqcurvefit, lsqnonlin) or a scalar (fmincon)
    if strcmp(this.fitopts.solverName, 'fmincon')
        out = sum(out.^2) ; 
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Laguerre polynomial L_{1/2}(x)
function l = Lhalf(x)
% Calculate Laguerre polynomial L_{1/2}(x)
% see Moments section of http://en.wikipedia.org/wiki/Rice_distribution

    l = exp(x/2) .* ( (1-x) .* besseli(0, -x/2) - x.*besseli(1, -x/2) ) ; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computes the log likelihood of the measurements given the model signals for the Rician noise model.
function loglikelihood = LogLikelihood(S, E, sigma)
% INPUTS: 
%       S: the measurements
%       E: computed signals from a model
%       sigma: the standard deviation of the Gaussian distributions underlying the Rician distribution.
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Aadapted from Alexander DC. Magn Reson Med. 2008; 60:439-448. 
% -------------------------------------------------------------------------------------------------------------------------
% Author: Junzhong Xu, 20160702
% -------------------------------------------------------------------------------------------------------------------------

    % check inputs
    if ~isscalar(sigma) , error('%s: The input sigma must be a scalar',mfilename) ; end
    
    % calculation
    loglikelihood = log(S) - log(sigma.^2) + logbesseli0(S.*E./(sigma.^2)) - (E.^2 + S.^2)./(2*sigma.^2) ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computes log(besseli(0,x)) robustly.  
function out = logbesseli0(x)
% Computing it directly causes % numerical problems at high x, but the function has asymptotic linear behaviour, which we approximate here for high x.
% For very large arguments to besseli, we approximate log besseli using a linear model of the asymptotic behaviour.
% The linear parameters come from this command: app=regress(log(besseli(0,500:700))',[ones(201,1) (500:700)']);
%  the calculation uses two segments:   lb0 = log(besseli(0,x)).*(x<700)+(x - log(2*pi*x)/2).*(x>=700) ; 
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Aadapted from the NODDI code (Author: Daniel C Alexander (d.alexander@ucl.ac.uk))
%       2. Revised by Junzhong Xu, 20160702
% -------------------------------------------------------------------------------------------------------------------------

    lb0 = zeros(length(x(:)), 1);
    indExact = find(x<700);
    indApprox = find(x>=700);
    % exact solution
    lb0(indExact) = log(besseli(0, x(indExact)));
    % This is a more standard approximation.  For large x, I_0(x) -> exp(x)/sqrt(2 pi x).
    lb0(indApprox) = x(indApprox) - log(2*pi*x(indApprox))/2 ;
    % reshape output
    out = reshape(lb0, size(x)) ; 

end
