function output = Signal_vin_d_Dex(parms, model)
%% IMPULSED with Din=const, betaex=const with model parms = [vin, d, Dex]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {vin, d, Dex} and will perform meshgrid(vin,d, Dex) to get all combinations of vin, d, and Dex
%               [option#2]: a 3xN numeric array. Necessary for data fitting. !!!! Note: rows are varying pulse parameters; columns are varying structure parms !!!
%       model: signal model
% OUTPUTS
%       out:  calculated dMRI signal
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Rows are varying pulse parameters; columns are varying structual parms
%       2. Make sure pulse (column vector) .* parms (row vector) to match dimensions, e.g., b(column) .* Dex (row)
% -------------------------------------------------------------------------------------------------------------------------

    %% check inputs
    if ~strcmp(model.structure.modelName, 'impulsed_vin_d_Dex'), error('%s: The input structure.ModelName is not ''impulsed_vin_d_Dex''',mfilename) ; end
    if ~isscalar(model.structure.Din), error('%s: model.structure.Din should a constant',mfilename); end
    if ~isscalar(model.structure.betaex), error('%s: model.structure.betaex should a constant',mfilename); end

    if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
        error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
    end
    if iscell(parms)            % input cell array
        for n = 1:length(parms)
            if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
        end
        if length(parms) ~= 3, error('%s: The input parms should take three (vin, d, Dex) parameters',mfilename) ; end
        % set up parms dimensions
        vin = parms{1}(:) ; d = parms{2}(:) ; Dex = parms{3}(:) ; 
        [vin, d, Dex] = meshgrid(vin, d, Dex) ; 
        parms = [vin(:)'; d(:)'; Dex(:)'] ; 
    elseif isnumeric(parms)    % input matrix
        if size(parms,1) ~=3, error('%s: The input parms should be a 2 x N numetric array', mfilename) ; end
    else
        error('%s: The input parms should be either a cell array or a 2xN numetric array',mfilename) ; 
    end
    
    % set Din the same dimension
    Din = repmat(model.structure.Din, [1, size(parms,2)]) ;
    
    %% calculation
    output = parms(1,:).*mig.Physics.RestrictedDWISignal([parms(2,:); Din], model) ...
        + (1-parms(1,:)).*exp(-model.pulse.b.*(parms(3,:) + model.structure.betaex'*model.pulse.f)) ; 
   
end
    