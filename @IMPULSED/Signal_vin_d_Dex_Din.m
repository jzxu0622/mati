function output = Signal_vin_d_Dex_Din(parms, model)
%% IMPULSED with betaex=const with model parms = [vin, d, Dex, Din]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {vin, d, Dex, Din} and will perform ndgrid(vin,d, Dex, Din) to get all combinations of vin, d, Dex, and Din
%               [option#2]: a 4xN numeric array. Necessary for data fitting. !!!! Note: rows are varying pulse parameters; columns are varying structure parms !!!
%       model: signal model
% OUTPUTS
%       out:  calculated dMRI signal
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Rows are varying pulse parameters; columns are varying structual parms
%       2. Make sure pulse (column vector) .* parms (row vector) to match dimensions, e.g., b(column) .* Dex (row)
% -------------------------------------------------------------------------------------------------------------------------


    %% check inputs
    if ~strcmp(model.structure.modelName, 'impulsed_vin_d_Dex_Din'), error('%s: The input structure.ModelName is not ''impulsed_vin_d_Dex_Din''',mfilename) ; end
    if ~isscalar(model.structure.betaex), error('%s: model.structure.betaex should a constant',mfilename); end

    if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
        error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
    end
    if iscell(parms)            % input cell array
        for n = 1:length(parms)
            if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
        end
        if length(parms) ~= 4, error('%s: The input parms should take three (vin, d, Dex, Din) parameters',mfilename) ; end
        % set up parms dimensions
        vin = parms{1}(:) ; d = parms{2}(:) ; Dex = parms{3}(:) ; Din = parms{4}(:) ; 
        [vin, d, Dex, Din] = ndgrid(vin, d, Dex, Din) ; 
        parms = [vin(:)'; d(:)'; Dex(:)'; Din(:)'] ; 
    elseif isnumeric(parms)    % input matrix
        if size(parms,1) ~= 4, error('%s: The input parms should be a 4 x N numetric array', mfilename) ; end
    else
        error('%s: The input parms should be either a cell array or a 4xN numetric array',mfilename) ; 
    end
    
    
    %% calculation
    output = parms(1,:).*mig.Physics.RestrictedDWISignal([parms(2,:); parms(4,:)], model) ...
        + (1-parms(1,:)).*exp(-model.pulse.b.*(parms(3,:) + model.pulse.f.*model.structure.betaex)) ; 
    
end
    