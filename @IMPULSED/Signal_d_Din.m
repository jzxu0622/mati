function out = Signal_d_Din(parms, model)
%% Calculate 1compartment restricted dMRI with model parms = [vin, d, Dex, Din]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {d, Din} and will perform meshgrid(d, Din) to get all combinations of d and Din
%               [option#2]: a 2xN numeric array. Necessary for data fitting. !!!! Note: rows are varying pulse parameters; columns are varying structure parms !!!
%       model: signal model
% OUTPUTS
%       out:  calculated dMRI signal
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. Rows are varying pulse parameters; columns are varying structual parms
%       2. Make sure pulse (column vector) .* parms (row vector) to match dimensions, e.g., b(column) .* Dex (row)
%       3. The structure geometry must be a member of {'plane','cylinder','sphere'}. The hollowSphere geometry is in another function mati.IMPULSED.Sigal_d_Din_hollow()
% -------------------------------------------------------------------------------------------------------------------------


    %% check inputs
    if ~strcmp(model.structure.modelName, '1compt'), error('%s: The input structure.ModelName is not ''impulsed_vin_d_Dex''',mfilename) ; end
    if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
        error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
    end
    if iscell(parms)            % input cell array
        for n = 1:length(parms)
            if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
        end
        if length(parms)~= 2, error('%s: The input parms should take two (d, Din) parameters',mfilename) ; end
        % set up parms dimensions
        d = parms{1}(:) ; D = parms{2}(:) ; 
        [d, D] = meshgrid(d, D) ; parms = [d(:)'; D(:)'] ; 
    elseif isnumeric(parms)    % input matrix
        if size(parms,1) ~=2, error('%s: The input parms should be a 2 x N numetric array', mfilename) ; end
    else
        error('%s: The input parms should be either a cell array or a 2xN numetric array',mfilename) ; 
    end
    
    %% calculate signals
    out = mati.Physics.RestrictedDWISignal(parms, model) ;
    
end