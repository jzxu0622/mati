classdef MRICytometry < mig.SignalModel
%% class definition of MRI-cytometry
% MRI-cytometry can take different structure models
% -------------------------------------------------------------------------------------------------------------------------
% Model choices: 
%           '1compt':       a single restricted geometry with a size d distribution and an intracellular diffusion coefficient Din distribution
%           'mricytometry_vin_d_Dex:   two-compartment with vin (intracellular proton fraction), d (intracelular size), and extracelluar diffusion coefficient Dex
%           'mricytometry_vin_d_Dex_Din:    two-compartment with vin, d, Dex, Din
%           'mricytometry_vin_d_Dex_Din_betaex:    two-compartment with vin, d, Dex, Din, betaex (extracellular diffusion linear dispersion rate)
% -------------------------------------------------------------------------------------------------------------------------
% NOTES: 
%       1. 1comptHollow is under construction. 
% -------------------------------------------------------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------------------------------------------------------

    %% properties
    properties (SetAccess = private)
        structure
        defaultStructure
    end
    
    %% static methods
    methods (Static)
        out = CreateDictionary(parms,this) ; 
        [x, resnorm] = FitFcn(this, signal, sigma, NumStarts) ;     % MRICytometry needs special fitting
    end

    %% other methods
    methods
        %%%%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%%%%
        % MRI-cytometry should take two inputs: a struct (structure model), an object (DiffusionPulseSequence), and an optional struct (fitopts)
        function this = MRICytometry(varargin)
            switch nargin
                case 0 , return ; 
                case 1 , error('%s: MRICytometry takes an object (DiffusionPulseSequence), a struct(structure),and [optional] struct (defaultFitopts)',mfilename) ; 
                case 2 , this.pulse = varargin{1} ; inputStructure = varargin{2} ; inputDefaultFitopts = [] ; 
                case 3 , this.pulse = varargin{1} ; inputStructure = varargin{2} ; inputDefaultFitopts = varargin{3} ; 
                otherwise , error('%s: MRICytometry takes an object (DiffusionPulseSequence), a struct(structure),and [optional] struct (defaultFitopts)',mfilename) ; 
            end

            %------------------------ setting defaultStructure --------------------------
            this.defaultStructure.modelName = [] ; this.defaultStructure.geometry = 'sphere' ; this.defaultStructure.Ndim = 3 ; 
            this.defaultStructure.Din = 1.56 ; this.defaultStructure.betaex = 0 ; this.defaultStructure.Dfree = 3.07 ; 
            this.parmsLabels = ["dvw","Dinvw","Dex","betaex","d","Din","dCal","DinCal","","","matrixInvw","matrixEx","matrixIn","matrixInCal"] ; this.NumParms = length(this.parmsLabels) ; 

            %------------------------ setting defaultFitopts -----------------------------
            this.defaultFitopts.dmin = 0.5 ; this.defaultFitopts.Dd = 0.5 ; this.defaultFitopts.dmax = 25 ; 
            this.defaultFitopts.Dinmin = 0.2 ; this.defaultFitopts.DDin = 0.2 ; this.defaultFitopts.Dinmax = this.defaultStructure.Dfree*1.1 ; 
            this.defaultFitopts.Dexmin = 0.2 ; this.defaultFitopts.DDex = 0.2 ; this.defaultFitopts.Dexmax = this.defaultStructure.Dfree*1.1 ; 
            this.defaultFitopts.betaexmin = 0 ; this.defaultFitopts.Dbetaex = 0.5 ; this.defaultFitopts.betaexmax = 10 ; 
            this.defaultFitopts.fittingMethod = 'special' ;          % {'normal','special'}. MRICytometry needs special fitting
            this.defaultFitopts.FitFcn = @this.FitFcn ; 
            this.defaultFitopts.penalty = 1e-2 ; this.defaultFitopts.penalty2 = 1e-4 ; 
            this.defaultFitopts.noiseModel = 'none' ;       % {'none','simple','rician','meanrician'}
            this.defaultFitopts.options =optimset('Display','off') ; 
            
            % ------------------ update structure using inputs -----------------------
            this.structure = this.defaultStructure ; 
            if ~isempty(inputStructure)
                names = fieldnames(inputStructure) ; 
                for n=1:length(names)
                    if ~isfield(this.structure, names{n}) , error('%s: The input structure has a field that is not a field of mig.MRICytometry.structure',mfilename) ; end
                    this.structure.(names{n}) = inputStructure.(names{n}) ; 
                end
            end
            % ------------------ update defaultFitopts using inputs, if any -----------------------
            if ~isempty(inputDefaultFitopts)
                names = fieldnames(inputDefaultFitopts) ; 
                for n=1:length(names)
                    if ~isfield(this.defaultFitopts, names{n}) , error('%s: The input defaultFitopts has a field that is not a field of mig.MRICytometry.defaultFitopts',mfilename) ; end
                    this.defaultFitopts.(names{n}) = inputDefaultFitopts.(names{n}) ; 
                end
            end
            
            % -------------------- create dictionary ---------------------------
            this.defaultFitopts = mig.MRICytometry.CreateDictionary(this) ; 
            
            % ------------------- other defaltFitopts -----------------------
            this.defaultFitopts.ADCfit_bmin = 0.19 ;                % minimum b for ADC fitting
            this.defaultFitopts.ADCfit_bmax = 1.05 ;             % maximum b for ADC fitting
            this.defaultFitopts.NumStarts = 50 ; 
            % flags for many options
            this.defaultFitopts.flag.denoise = 'n' ; 
            this.defaultFitopts.flag.degibbs = 'n' ; 
            this.defaultFitopts.flag.deivim = 'y' ; 
            this.defaultFitopts.flag.multistart = 'n' ; 
            this.defaultFitopts.flag.parfor = 'n' ; 
            
        end % End of Constructor
        
    end % End of Methods
    
end % End of class
