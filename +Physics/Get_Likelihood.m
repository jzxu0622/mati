function likelihood = get_Likelihood(S, E, sigma)
%% Computes the log likelihood of the measurements given the model signals for the Rician noise model.
%
% loglikelihood = MRI_LogLikelihood(S, E, sigma) returns the likelihood of measuring meas given the signals and the noise standard deviation sig. 
% Inputs: 
% S: the measurements
% E: computed signals from a model
% sigma: the standard deviation of the Gaussian distributions underlying the Rician distribution.
%
% ------------------------------------------------------
% Originally written by Daniel C Alexander (d.alexander@ucl.ac.uk)
% 
% Aadapted from NODDI matlab code
% jzxu 20160702
% There was a bug in original code: the first term should be log(S), not log(E) 
% ------------------------------------------------------

i_nonzero=find(S~=0);
S=S(i_nonzero);
E=E(i_nonzero);

if length(S)~=length(E(:))
    sprintf('the dimension of measured signals does not match the dimension of fitted signals...')
    likelihood=0;
    return;
%     S = expmatrix(E*0,S,length(size(E)));
end
S = double(S);
logliks = log(S) - log(sigma.^2) + logbesseli0(S.*E./(sigma.^2)) - (E.^2 + S.^2)./(2*sigma.^2) ;
loglikelihood = sum(logliks);
likelihood = exp(loglikelihood);
%likelihood = likelihood/sum(likelihood(:)); 
end

function lb0 = logbesseli0(x)
%% Computes log(besseli(0,x)) robustly.  
% Computing it directly causes % numerical problems at high x, but the function has asymptotic linear behaviour, which we approximate here for high x.
%
% author: Daniel C Alexander (d.alexander@ucl.ac.uk)
%

% For very large arguments to besseli, we approximate log besseli using a linear model of the asymptotic behaviour.
% The linear parameters come from this command: 
% app=regress(log(besseli(0,500:700))',[ones(201,1) (500:700)']);
    x_size = size(x) ;     
    app = [-3.61178295877576 0.99916157999904];
    lb0 = zeros(length(x(:)), 1);
    exact = find(x<700);
    approx = find(x>=700);
    lb0(exact) = log(besseli(0, x(exact)));
    %lb0(approx) = x(approx)*app(2) + app(1);
    % This is a more standard approximation.  For large x, I_0(x) -> exp(x)/sqrt(2 pi x).
    lb0(approx) = x(approx) - log(2*pi*x(approx))/2 ;
    lb0 = reshape(lb0, x_size) ;     
%     lb0 = log(besseli(0,x)).*(x<700)+(x - log(2*pi*x)/2).*(x>=700) ; 

end