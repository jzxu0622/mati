%% Example of performing IMPULSED fitting
%
% -------------------------------------------------------------------------------------------------------------------------

%% preliminary
clear variables ; clear obj ; 

%% generate DiffusionPulseSequence object
pulse_tcos = mig.DiffusionPulseSequence(9,...
    'TE',               110,...
    'delta',            40,...
    'Delta',            45,...
    'b',                  [0.1,0.2,0.3, 0.4, 0.2,0.4,0.6,0.8,1.0], ...
    'n',                  [2,2,2, 2, 1,1,1,1,1],...
    'shape',         "tcos",...
    'gdir',             [0 0 1],...
    'trise',            0.9) ; 

pulse_tpgse = mig.DiffusionPulseSequence(9,...
    'TE',               110,...
    'delta',            12, ...
    'Delta',            74, ...
    'b',                  [0.2:0.2:1.8], ...
    'shape',         "tpgse",...
    'gdir',             [0 0 1],...
    'trise',             0.9) ; 

pulse = mig.PulseSequence.cat(pulse_tcos, pulse_tpgse) ; 

% fitting range
pulse = pulse(pulse.G<80e-5) ; % 80mT/m 
pulse.disp(pulse)


%% Generate IMPULSED model object
% choose which model to use. All structure fields set here will override the default structure fields. 
nmodel = 3 ; 
switch nmodel
    case 1, structure.modelName = '1compt' ; structure.geometry = 'sphere' ; 
    case 2, structure.modelName = 'impulsed_vin_d_Dex' ; structure.Din = 2 ; structure.betaex = 0 ; structure.geometry = 'sphere';
    case 3, structure.modelName = 'impulsed_vin_d_Dex_Din' ; %structure.betaex = 0 ; structure.geometry = 'sphere';
    case 4, structure.modelName = 'impulsed_vin_d_Dex_Din_betaex' ; %structure.geometry = 'sphere';
end
impulsed = mig.IMPULSED(structure, pulse) ; 


%% [Optional] Test IMPULSED model
% generate my test structure. rows are #of unique parms, columns are different parms
switch nmodel
    case 1      % [d, Din]
        d = [10:15] ; Din = [1.56 3] ; 
        parms_sim = {d, Din};    [d,Din]=meshgrid(d,Din) ; d_sim = d(:)' ; Din_sim = Din(:)' ; 
    case 2      % [vin, d, Dex]
        vin = [0.6] ; d = [10:15] ; Dex = [1.56 3] ; 
        parms_sim = {vin, d, Dex};   [vin, d,Dex] = meshgrid(vin, d,Dex) ; d_sim = d(:)' ; 
    case 3      % [vin, d, Dex, Din]
        vin = [0.6] ; d = [10:15] ; Dex = [1.56 3] ; Din = [1.56] ; 
        parms_sim = {vin, d, Dex, Din};   [vin, d,Dex,Din] = ndgrid(vin, d,Dex,Din) ; d_sim = d(:)' ; 
%         [vin, d,Dex, Din] = ndgrid(vin, d,Dex,Din) ; parms_sim = [vin(:)'; d(:)'; Dex(:)'; Din(:)'] ;   d_sim = d(:)' ; 
    case 4      % [vin, d, Dex, Din, betaex]
        vin = [0.6] ; d = [8:2:16] ; Dex = [2] ; Din = [1.56] ; betaex = [5] ; 
        parms_sim = {vin, d, Dex, Din, betaex};   [vin, d,Dex,Din,betaex] = ndgrid(vin, d,Dex,Din,betaex) ; d_sim = d(:)' ; 
end

signal_sim = impulsed.FcnSignal(parms_sim, impulsed) ; 
signal_sim = mig.Physics.AddRicianNoise(signal_sim,0.025) ; 
[Npulse, Nparms] = size(signal_sim) ; 
data = mig.ImageData(reshape(signal_sim',[Nparms, 1 1 Npulse]),0.025) ; 


%% Fit signal
% create fit object
fitopts.solverName = 'lsqnonlin'; % {'lsqcurvefit' , 'lsqnonlin' , 'fmincon'}
fitopts.options = optimoptions(fitopts.solverName,'Display','off') ; 
fitopts.noiseModel = 'none' ; %{'none','simple';'rician'}
fitopts.flag.multistart = 'y' ; fitopts.flag.parfor = 'y' ; fitopts.flag.deivim = 'n' ; 
fitopts.NumStarts = 5 ; 

fitpars = mig.FitPars(impulsed, fitopts) ; 
warning off ; 
% fit model to data
fitout = fitpars.Fit(data) ; 


%% check fitted results
figure(1) ; clf ; hold on ; 
plot(d_sim, fitout.d, 'o') ; 
plot([0 20],[0 20],'r') ; box on ; 
xlabel('input d [\mum]') ; ylabel('fitted d [\mum]') ; xlim([0 20]) ; ylim([0 20]) ; 
legend('fits', 'identity', 'Location','Southeast') ; 


%% End