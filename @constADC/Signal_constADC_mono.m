function out = Signal_constADC_mono( parms, model )
%% model parms = [S0, ADC]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {vin, d, Dex, kin} and will
%               perform ndgrid(vin,d, Dex, kin) to get all combinations of vin,d, Dex, Din, and kin
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
if ~strcmp(model.structure.modelName, 'constADC_mono'), error('%s: The input structure.ModelName is not ''constADC_mono''',mfilename) ; end
if iscell(parms)            % input cell array
    for n = 1:length(parms)
        if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
    end
    if length(parms) ~= 2, error('%s: The input parms should take two (S0, ADC) parameters',mfilename) ; end
    % set up parms dimensions
    S0 = parms{1}(:) ; ADC = parms{2}(:) ;  
    [S0, ADC] = ndgrid(S0, ADC) ; 
    parms = [S0(:)'; ADC(:)'] ; 
elseif isnumeric(parms)    % input matrix
    if size(parms,1) ~= 2, error('%s: The input parms should be a 2 x N numetric array', mfilename) ; end
    S0 = parms(1,:);
    ADC = parms(2,:);    
else
    error('%s: The input parms should be either a cell array or a 2xN numetric array',mfilename) ;
end
%% calculation        
out=S0.*exp(-model.pulse.b.*ADC);
out(model.pulse.b<1e-4)=1;
end