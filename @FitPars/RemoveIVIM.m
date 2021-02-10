function out = RemoveIVIM(this, data)
%% mig.IMPULSED.RemoveIVIM to remove IVIM effects in diffusion weighted data
% INPUTS:
%       this: a mig.FitPars object
%       data: a mig.ImageData object
% OUTPUTS:
%       out: a mig.ImageData object
% -------------------------------------------------------------------------------------------------------------------------

    %% preliminary
    warning('off') ; 
    
    % check input 
    if ~isa(data,'mig.ImageData'), error('%s: the input data should be an ImageData object',mfilename) ; end
    
    % preparation 
    img = zeros([data.Nx, data.Ny, data.Nz, data.Nt]) ; 
    pulse = this.model.pulse ; 
    
    % options
    if this.fitopts.flag.parfor == 'y'
        poolobj = parpool ;  Ncore = poolobj.NumWorkers ; 
    else
        Ncore = 0 ; 
    end
        
        
    %% fitting
    % points
    [row, col, slice]=ind2sub(size(data.mask),find(data.mask>0)) ;     Npixel = length(row) ; 
    img1d = zeros([data.Nt, Npixel]) ;    
    
%     for npixel = 1:Npixel
    parfor (npixel = 1:Npixel, Ncore)
        nx = row(npixel) ; ny = col(npixel) ; nz = slice(npixel) ; 
        signal_raw = squeeze(data.img(nx,ny,nz,:)) ; 
        
        % preliminary 
        signal = zeros(size(signal_raw)) ; 
        tdiffs = unique(pulse.tdiff) ; Ntdiff = length(tdiffs) ; 
    
        for ntdiff = 1:Ntdiff
            ind_tdiff = (pulse.tdiff==tdiffs(ntdiff)) ; 
            ind_fit = (pulse.tdiff==tdiffs(ntdiff)) & (pulse.b > this.model.defaultFitopts.ADCfit_bmin) & (pulse.b < this.model.defaultFitopts.ADCfit_bmax) ; 
            switch sum(ind_fit)
                case {0,1}      % not enough points to perform ADC fitting so keep original signal
%                     warning on ; warning('%s: Not sufficient points for ADC fitting to remove IVIM') ; warning off ; 
                    ind_fit = ind_tdiff ; S0 = 1 ; 
                case {2,3}      % exponential ADC fitting, Gaussian Phase Approximation
                    p = polyfit(pulse.b(ind_fit), -log(signal_raw(ind_fit)),1) ; 
                    S0 = exp(-polyval(p,0)) ; 
                otherwise       % biexponential fitting to encounter for any non-Gaussian diffusion 
                     [x, ~,~] = MRI_fit_biexp(pulse.b(ind_fit), signal_raw(ind_fit), [1,0.2,0.5,1], [ 0 0 0 0], [2 4 4 1]) ; 
                     S0 = x(4) ; 
            end
            
            % removal IVIM
            signal(ind_tdiff) = signal_raw(ind_tdiff) / S0 ; 
        end
                   
        img1d(:, npixel) = signal ; 
        
    end % end of parfor loop
    
    %% output results
    for npixel = 1:Npixel
        img(row(npixel),col(npixel),slice(npixel),:) = reshape(img1d(:,npixel), [1 1 1 size(img1d,1)]) ;  
    end
    
    out = data.UpdateImg(img) ; 
    
end  % end of function
