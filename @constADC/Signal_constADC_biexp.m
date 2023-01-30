function out = Signal_constADC_biexp( parms, model )
%% Model including both diffusion and water exchange with all five parameters with model parms = [vin, d, Dex, Din, kin]
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
if ~strcmp(model.structure.modelName, 'constADC_biexp'), error('%s: The input structure.ModelName is not ''constADC_biexp''',mfilename) ; end

if iscell(parms)            % input cell array
    for n = 1:length(parms)
        if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
    end
    if length(parms) ~= 4, error('%s: The input parms should take four (S0, f, ADC1, ADC2) parameters',mfilename) ; end
    % set up parms dimensions
    S0 = parms{1}(:) ; f = parms{2}(:) ; ADC1 = parms{3}(:) ; ADC2 = parms{4}(:) ; 
    [S0, f, ADC1, ADC2] = ndgrid(S0, f, ADC1, ADC2) ; 
    parms = [S0(:)'; f(:)'; ADC1(:)'; ADC2(:)'] ; 
elseif isnumeric(parms)    % input matrix
    if size(parms,1) ~= 4, error('%s: The input parms should be a 4 x N numetric array', mfilename) ; end
    S0 = parms(1,:);
    f = parms(2,:);
    ADC1 = parms(3,:);    
    ADC2 = parms(4,:);
else
    error('%s: The input parms should be either a cell array or a 4xN numetric array',mfilename) ;
end
%% calculation  
out=s0.*(f.*exp(-model.pulse.b.*ADC1)+(1-f).*exp(-model.pulse.b.*ADC2));
out(model.pulse.b<1e-4)=1;
end