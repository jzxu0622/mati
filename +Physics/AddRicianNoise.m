function dataNoisy = AddRicianNoise(data, sigma)
%% Add Rician noise to data
% INPUTS: 
%       data:        imaging data w/o noise
%       sigma:      noise sigma. It should be empty (no noise) or a scalar (all data points have the same noise
% OUTPUTS: 
%       dataNoisy: data with Rician noise
% ---------------------------------------------------------------------------------
% revised by Junzhong Xu, Jan 29, 2020
% 
% ---------------------------------------------------------------------------------


    %% check inputs
    if isempty(sigma)
        warning('%s: The noise sigma is empty so no noise is added to signal',mfilename) ; 
        dataNoisy = data ; 
    end
    
    if ~isscalar(sigma)
        error('%s: The noise sigma should be either empty (no noise) or a scalar (all data points have the same noise)',mfilename) ; 
    end
    
    %% add noise to data
    Sx = sigma .* randn(size(data)) + data ; 
    Sy = sigma .* randn(size(data)) ; 
    dataNoisy = sqrt(Sx.^2 + Sy.^2) ; 

end
