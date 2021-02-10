function out = Signal_vin_d_Dex_Din_kin( parms, model )
%% Model including both diffusion and water exchange with all five parameters with model parms = [vin, d, Dex, Din, kin]
% INPUTS
%       parms: Input parameters. 
%               [option#1]: a cell array {vin, d, Dex, Din, kin} and will
%               perform ndgrid(vin,d, Dex, Din, kin) to get all combinations of vin,d, Dex, Din, and kin
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
if ~strcmp(model.structure.modelName, 'exchange_vin_d_Dex_Din_kin'), error('%s: The input structure.ModelName is not ''exchange_vin_d_Dex_Din_kin''',mfilename) ; end
if ~ismember(model.structure.geometry, {'plane','cylinder','sphere'})
    error('%s: The input structure.geometry must be plane, cylinder, or sphere geometry',mfilename) ; 
end
if iscell(parms)            % input cell array
    for n = 1:length(parms)
        if ~isnumeric(parms{n}), error('%s: The input parms\{%d\} should be numeric',mfilename,n); end     
    end
    if length(parms) ~= 5, error('%s: The input parms should take five (vin, d, Dex, Din, kin) parameters',mfilename) ; end
    % set up parms dimensions
    vin = parms{1}(:) ; d = parms{2}(:) ; Dex = parms{3}(:) ; Din = parms{4}(:) ; kin = parms{5}(:) ; 
    [vin, d, Dex, Din, kin] = ndgrid(vin, d, Dex, Din, kin) ; 
    parms = [vin(:)'; d(:)'; Dex(:)'; Din(:)'; kin(:)'] ; 
elseif isnumeric(parms)    % input matrix
    if size(parms,1) ~= 5, error('%s: The input parms should be a 5 x N numetric array', mfilename) ; end
    vin = parms(1,:);
    d = parms(2,:);
    Dex = parms(3,:);
    Din = parms(4,:);
    kin = parms(5,:);
else
    error('%s: The input parms should be either a cell array or a 5xN numetric array',mfilename) ;
end
%% calculation    

kex=vin.*kin./(1-vin);
kex(isinf(kex)|isnan(kex))=0;
mask_PGSE=(model.pulse.n==0);
mask_OGSE=(model.pulse.n>0);

imp_structure=model.structure;
imp_structure.modelName = 'impulsed_vin_d_Dex_Din';
imp_model = mig.IMPULSED(imp_structure, model.pulse) ; 
E = RestrictedDWISignal([parms(2,:); parms(4,:)], imp_model);
ADC = -log(E)./model.pulse.b;
% PGSE
tg1=(model.pulse.TE-model.pulse.tdiff)/2;
tg2=tg1+model.pulse.tdiff;
tg3=model.pulse.TE;

[ M_in_tg1, M_ex_tg1 ] = mig.Physics.generalSolution_blochEQ( kin,kex,kin,kex,vin,1-vin,tg1 );
Ain=kin+model.pulse.b./model.pulse.tdiff.*ADC;
Aex=kex+model.pulse.b./model.pulse.tdiff.*Dex;
[ M_in_tg2, M_ex_tg2 ] = mig.Physics.generalSolution_blochEQ( Ain,Aex,kin,kex,M_in_tg1,M_ex_tg1,tg2-tg1 );
[ M_in_tg3, M_ex_tg3 ] = mig.Physics.generalSolution_blochEQ( kin,kex,kin,kex,M_in_tg2,M_ex_tg2,tg3-tg2 );
signal_PGSE=(M_in_tg3+M_ex_tg3).*mask_PGSE;
        
% OGSE
tg1=model.pulse.TE/2-model.pulse.delta/4-model.pulse.Delta/2;
Ain=kin;
Aex=kex;
[ M_in_tg1, M_ex_tg1 ] = mig.Physics.generalSolution_blochEQ( Ain,Aex,kin,kex,vin,1-vin,tg1 );
Ain=kin+model.pulse.b./model.pulse.tdiff./model.pulse.n/4.*ADC;
Aex=kex+model.pulse.b./model.pulse.tdiff./model.pulse.n/4.*Dex;
[ M_in_tg2, M_ex_tg2 ] = mig.Physics.generalSolution_blochEQ( Ain,Aex,kin,kex,M_in_tg1,M_ex_tg1,model.pulse.delta/2);
[ M_in_tg3, M_ex_tg3 ] = mig.Physics.generalSolution_blochEQ( kin,kex,kin,kex,M_in_tg2,M_ex_tg2,model.pulse.Delta-model.pulse.delta/2);
[ M_in_tg4, M_ex_tg4 ] = mig.Physics.generalSolution_blochEQ( Ain,Aex,kin,kex,M_in_tg3,M_ex_tg3,model.pulse.delta/2);
[ M_in_tg5, M_ex_tg5 ] = mig.Physics.generalSolution_blochEQ( kin,kex,kin,kex,M_in_tg4,M_ex_tg4,tg1);
signal_OGSE=(M_in_tg5+M_ex_tg5).*mask_OGSE;
signal_OGSE(isnan(signal_OGSE)|isinf(signal_OGSE))=0;
signal_PGSE(isnan(signal_PGSE)|isinf(signal_PGSE))=0;
out=signal_PGSE+signal_OGSE;

end