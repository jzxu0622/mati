%% Example of performing MRICytometry fitting
%
% -------------------------------------------------------------------------------------------------------------------------

%% preliminary
clear variables ; clear obj ; 

%% generate DiffusionPulseSequence object
pulse_tcos = mati.DiffusionPulseSequence(9,...
    'TE',               110,...
    'delta',            40,...
    'Delta',            45,...
    'b',                  [0.1,0.2,0.3, 0.4, 0.2,0.4,0.6,0.8,1.0], ...
    'n',                  [2,2,2, 2, 1,1,1,1,1],...
    'shape',         "tcos",...
    'gdir',             [0 0 1],...
    'trise',            0.9) ; 

pulse_tpgse = mati.DiffusionPulseSequence(9,...
    'TE',               110,...
    'delta',            12, ...
    'Delta',            74, ...
    'b',                  [0.2:0.2:1.8], ...
    'shape',         "tpgse",...
    'gdir',             [0 0 1],...
    'trise',             0.9) ; 

pulse = mati.PulseSequence.cat(pulse_tcos, pulse_tpgse) ; 

% fitting range
pulse = pulse(pulse.G<80e-5) ; % 80mT/m 
pulse.disp(pulse)

%% Generate MRICytometry model object
% choose which model to use. All structure fields set here will override the default structure fields. 
nmodel = 7 ; 
switch nmodel
    case 1, structure.modelName = 'dNone_DinDist_DexNone' ; structure.Dfree = 3.07 ;  
    case 2, structure.modelName = 'dDist_DinFixed_DexNone' ;  structure.Din = 1.56 ;  fitopts.Dd=1 ; 
    case 3, structure.modelName = 'dDist_DinDist_DexNone' ;  fitopts.Dd=1 ; structure.Dfree = 3.07 ;  
    case 4, structure.modelName = 'dDist_DinFixed_DexDist' ;     structure.Din = 1.56 ; 
    case 5, structure.modelName = 'dDist_DinDist_DexDist' ;     
    case 6, structure.modelName = 'dDist_DinFixed_DexDisper' ;     
    case 7, structure.modelName = 'dDist_DinDist_DexDisper' ;     
    otherwise , error('The nmodel is not reconized') ; 
end

fitopts.fittingMethod = 'special' ; 
mricytometry = mati.MRICytometry(pulse, structure, fitopts) ; 


%% load data here. [Optional] Test MRICytometry model using simulated data
% user input parameters for simulations
nDist = 1 ; 
% set d distribution
switch nDist
    case 1
        sim.distribution = 'gaussian' ; 
        sim.dcen = 16 ; sim.dsigma = 3 ; 
    case 2
        sim.distribution = 'bi-modal' ; 
        sim.dcen1 = 8 ; sim.dsigma1 = 2 ; sim.frac1 = 0.7 ; 
        sim.dcen2 = 16 ; sim.dsigma2 = 2 ;         
    case 3
        sim.distribution = 'gamma' ; 
        sim.dalpha = 4 ; sim.dbeta = 2 ; 
    otherwise
end
% set other parms
sim.Ndim = 3 ; % 1d, 2d, or 3d
sim.Dincen = 1.58 ; sim.Dinsigma = 0.5 ; sim.Dexcen = 2 ; sim.Dexsigma = 0.5 ;  sim.betaexcen = 2; sim.betaexsigma = 0.25 ; 
sim.vin = 0.7; sim.vex = 0.3 ; sim.vfree = 1 - sim.vin - sim.vex ; 
sim.flag.DinDist = 'y' ; 
% calculate simulated parameters
sim = FcnSetSimPars(sim) ; 
% calculate simulated signal
signal_sim = FcnSimulateSignal(sim, mricytometry) ; 


%%  create ImageData object
noiseSigma = 1e-5 ; 
signal_sim = mati.Physics.AddRicianNoise(signal_sim, noiseSigma) ; 
[Npulse, Nparms] = size(signal_sim) ; 
data = mati.ImageData(reshape(signal_sim',[Nparms, 1 1 Npulse])) ; 


%% Fit signal
% create fit object
fitopts.noiseModel = 'none' ; %{'none','simple';'rician'}
fitopts.flag.multistart = 'y' ; fitopts.flag.parfor = 'n' ; fitopts.flag.deivim = 'n' ; 

fitpars = mati.FitPars(mricytometry, fitopts) ; 
warning off ; 
% fit model to data
fitout = fitpars.Fit(data) ; 

%% show results
FcnPlotDiagram(sim, fitout.parms{1,1,1}, fitpars.fitopts) ; 

%% End