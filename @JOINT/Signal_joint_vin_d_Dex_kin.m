function out = Signal_joint_vin_d_Dex_kin( parms, model )
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
if ~strcmp(model.structure.modelName, 'joint_vin_d_Dex_kin'), error('%s: The input structure.ModelName is not ''joint_vin_d_Dex_kin''',mfilename) ; end
if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
    error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
end
if iscell(parms)            % input cell array
    for n = 1:length(parms)
        if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
    end
    if length(parms) ~= 4, error('%s: The input parms should take four (vin, d, Dex, kin) parameters',mfilename) ; end
    % set up parms dimensions
    vin = parms{1}(:) ; d = parms{2}(:) ; Dex = parms{3}(:) ; kin = parms{4}(:) ; 
    [vin, d, Dex, kin] = ndgrid(vin, d, Dex, kin) ; 
    parms = [vin(:)'; d(:)'; Dex(:)'; kin(:)'] ; 
elseif isnumeric(parms)    % input matrix
    if size(parms,1) ~= 4, error('%s: The input parms should be a 4 x N numetric array', mfilename) ; end
    vin = parms(1,:);
    d = parms(2,:);
    Dex = parms(3,:);    
    kin = parms(4,:);
else
    error('%s: The input parms should be either a cell array or a 4xN numetric array',mfilename) ;
end
%% calculation  
% set Din the same dimension
Din = repmat(model.structure.Din, [1, size(parms,2)]) ;
kex=vin.*kin./(1-vin);
kex(isinf(kex)|isnan(kex))=0;
mask_PGSE=(model.pulse.n==0);
mask_OGSE=(model.pulse.n>0);

imp_structure=model.structure;
imp_structure.modelName = 'impulsed_vin_d_Dex';
imp_model = mati.IMPULSED(imp_structure, model.pulse) ; 
E = mati.Physics.RestrictedDWISignal([d; Din], imp_model);
ADC = -log(E)./model.pulse.b;
% PGSE
tg1=(model.pulse.TE-model.pulse.tdiff)/2;
tg2=tg1+model.pulse.tdiff;
tg3=model.pulse.TE;

[ M_in_tg1, M_ex_tg1 ] = mati.JOINT.generalSolution_blochEQ( kin,kex,kin,kex,vin,1-vin,tg1 );
Ain=kin+model.pulse.b./model.pulse.tdiff.*ADC;
Aex=kex+model.pulse.b./model.pulse.tdiff.*Dex;
[ M_in_tg2, M_ex_tg2 ] = mati.JOINT.generalSolution_blochEQ( Ain,Aex,kin,kex,M_in_tg1,M_ex_tg1,tg2-tg1 );
[ M_in_tg3, M_ex_tg3 ] = mati.JOINT.generalSolution_blochEQ( kin,kex,kin,kex,M_in_tg2,M_ex_tg2,tg3-tg2 );
signal_PGSE=(M_in_tg3+M_ex_tg3).*mask_PGSE;
        
% OGSE
signal_OGSE=vin.*E + (1-vin).*exp(-model.pulse.b.*(Dex + model.pulse.f.*model.structure.betaex));
signal_OGSE=signal_OGSE.*mask_OGSE;
signal_OGSE(isnan(signal_OGSE)|isinf(signal_OGSE))=0;
signal_PGSE(isnan(signal_PGSE)|isinf(signal_PGSE))=0;
out=signal_PGSE+signal_OGSE;
out(model.pulse.b<1e-4)=1;
end