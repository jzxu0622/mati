function out = Signal_d_Din_hollow(parms, model)
%% Calculate 1compt DWI signal restricted in a hollow sphere. Model parms=[din, dout, Din]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {din, dout, Din} and will perform ndgrid(din,dout, Din) to get all combinations of din, dout, and Din
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
    if ~strcmp(model.structure.modelName, '1comptHollow'), error('%s: The input structure.ModelName is not ''1comptHollow''',mfilename) ; end
    if ~ismember(model.structure.geometry, {'hollowSphere','sphericalShell'})
        error('%s: The input structure.geometry must be hollowSphere or sphericalShell',mfilename) ; 
    end
    if iscell(parms)            % input cell array
        for n = 1:length(parms)
            if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
        end
        if length(parms)~= 3, error('%s: The input parms should take two (d, Din) parameters',mfilename) ; end
        % set up parms dimensions
        din = parms{1}(:) ; dout = parms{2}() ; D = parms{3}(:) ; 
        [din, dout, D] = meshgrid(din, dout, D) ; parms = [din(:)'; dout(:)' ; D(:)'] ; 
    elseif isnumeric(parms)    % input matrix
        if size(parms,1) ~= 3, error('%s: The input parms should be a 3 x N numetric array', mfilename) ; end
    else
        error('%s: The input parms should be either a cell array or a 3xN numetric array',mfilename) ; 
    end
    
    %% calculation
    out = mig.Physics.RestrictedDWISignal(parms, model) ;
    
end