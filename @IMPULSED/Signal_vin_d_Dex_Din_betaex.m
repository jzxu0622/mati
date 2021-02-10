function out = Signal_vin_d_Dex_Din_betaex(parms, model)
%% IMPULSED with all five parameters with model parms = [vin, d, Dex, Din, betaex]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {vin, d, Dex, Din, betaex} and will perform ndgrid(vin,d, Dex, Din, betaex) to get all combinations of vin,d, Dex, Din, and betaex
%               [option#2]: a 5xN numeric array. Necessary for data fitting. !!!! Note: rows are varying pulse parameters; columns are varying structure parms !!!
%       model: signal model
% OUTPUTS
%       out:  calculated dMRI signal
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Rows are varying pulse parameters; columns are varying structual parms
%       2. Make sure pulse (column vector) .* parms (row vector) to match dimensions, e.g., b(column) .* Dex (row)
% -------------------------------------------------------------------------------------------------------------------------

    %% check inputs
    if ~strcmp(model.structure.modelName, 'impulsed_vin_d_Dex_Din_betaex'), error('%s: The input structure.ModelName is not ''impulsed_vin_d_Dex_Din_betaex''',mfilename) ; end
    if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
        error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
    end
    if iscell(parms)            % input cell array
        for n = 1:length(parms)
            if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
        end
        if length(parms) ~= 5, error('%s: The input parms should take three (vin, d, Dex, Din, betaex) parameters',mfilename) ; end
        % set up parms dimensions
        vin = parms{1}(:) ; d = parms{2}(:) ; Dex = parms{3}(:) ; Din = parms{4}(:) ; betaex = parms{5}(:) ; 
        [vin, d, Dex, Din, betaex] = ndgrid(vin, d, Dex, Din, betaex) ; 
        parms = [vin(:)'; d(:)'; Dex(:)'; Din(:)'; betaex(:)'] ; 
    elseif isnumeric(parms)    % input matrix
        if size(parms,1) ~= 5, error('%s: The input parms should be a 5 x N numetric array', mfilename) ; end
    else
        error('%s: The input parms should be either a cell array or a 5xN numetric array',mfilename) ; 
    end
    
    %% calculation
    out = parms(1,:).*mati.Physics.RestrictedDWISignal([parms(2,:); parms(4,:)], model) ...
        + (1-parms(1,:)).*exp(-model.pulse.b.*(parms(3,:) + model.pulse.f.*parms(5,:)));
    
end
    