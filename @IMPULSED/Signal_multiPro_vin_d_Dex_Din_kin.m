function out = Signal_multiPro_vin_d_Dex_Din_kin( parms, model )
%% Multi-propagator Model including both diffusion and water exchange with all five parameters with model parms = [vin, d, Dex, Din, kin]
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
if ~strcmp(model.structure.modelName, 'multiPro_vin_d_Dex_Din_kin'), error('%s: The input structure.ModelName is not ''multiPro_vin_d_Dex_Din_kin''',mfilename) ; end
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
    vin = parms(1,1);
    d = parms(2,1);
    Dex = parms(3,1);
    Din = parms(4,1);
    kin = parms(5,1);
else
    error('%s: The input parms should be either a cell array or a 5xN numetric array',mfilename) ;
end
if ~isfield(model,'propagatorSize')
    matrixSize=model.defaultPropagatorSize;  % matrix size for the propagators
else
    matrixSize=model.propagatorSize;
end
%% calculation    
kex=kin.*vin./(1-vin);
kex(isnan(kex)|isinf(kex))=0;
timeStep=model.pulseArray.delta_t;
out = zeros(length(vin(:)),model.pulse.Nacq);
for i_n=1:length(vin(:))
    beta=[vin(i_n) d(i_n) Dex(i_n) Din(i_n) kin(i_n) kex(i_n)];
    R=mig.Physics.calMatrixR(timeStep,matrixSize,model.structure,beta);
    for i_b=1:model.pulse.Nacq
        [i_n i_b]
        n_pulse=model.pulseArray.n(i_b)-1;
        if model.pulse.b(i_b)==0 
            out(i_n,i_b)=1;
            continue;
        end
        Dex_t=model.pulse.b(i_b)*Dex(i_n)/(model.pulse.delta(i_b)+model.pulse.Delta(i_b));  % if we include a diffusion tim -dependency for Dex, we need to modify the expression for Dex here.
        s_ex=1-vin(i_n);
        A_0=mig.Physics.calMatrixA(model.pulse.gamma/2/pi*model.pulse.G(i_b)*timeStep/max(model.pulseArray.delta_q(i_b,:)),matrixSize,model.structure,beta); % A with the smallest q
        A_max=mig.Physics.calMatrixA(model.pulse.gamma/2/pi*model.pulse.G(i_b)*timeStep,matrixSize,model.structure,beta);
        
        for i=1:n_pulse
            if i==1
                S=mig.Physics.calMatrixS(model.pulseArray.q(i_b,i)*model.pulse.gamma/2/pi, matrixSize,model.structure,beta);       
                s_s=vin(i_n)*S;
                continue;
            end
            if i==n_pulse
                S=mig.Physics.calMatrixS(abs(model.pulseArray.q(i_b,i))*model.pulse.gamma/2/pi, matrixSize,model.structure,beta);
                S_nag=ctranspose(S);
                s_s=s_s*R*S_nag;
                continue;
            end            
            s_in=abs(sum(s_s(:)));     
            if model.pulseArray.delta_q(i_b,i)>0
                if (model.pulseArray.delta_q(i_b,i)-max(model.pulseArray.delta_q(i_b,:)))<1e-2
                    A=A_max;
                else
                    A=1;
                    for ii=1:model.pulseArray.delta_q(i_b,i)
                        A=A*A_0;
                    end            
                end
                s_s=s_s*R*A;
            elseif model.pulseArray.delta_q(i_b,i)<0
                if (abs(model.pulseArray.delta_q(i_b,i))-max(model.pulseArray.delta_q(i_b,:)))<1e-2
                    A=A_max;
                else
                    A=1;
                    for ii=1:abs(model.pulseArray.delta_q(i_b,i))
                        A=A*A_0;
                    end            
                end
                A_nag=ctranspose(A);
                s_s=s_s*R*A_nag;
            else
                s_s=s_s*R;
            end
            Dinstant=-log(abs(sum(s_s(:)))/s_in)/timeStep;
            [s_in_t,s_ex]=mig.Physics.generalSolution_blochEQ( Dinstant+kin(i_n),Dex_t+kex(i_n),kin(i_n),kex(i_n),s_in,s_ex,timeStep); 
            s_s=s_s/abs(sum(s_s(:)))*s_in_t;
        end       
        out(i_n,i_b)=abs(sum(s_s(:)))+s_ex;
    end    
end
end