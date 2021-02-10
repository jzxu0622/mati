classdef IMPULSED < mig.SignalModel
%% class definition of IMPULSED  
% IMPULSED can take different structure models
% Model choices: 
%           '1compt':       a single restricted geometry with a size d and an intracellular diffusion coefficient Din
%           '1comptHollow':     a single restricted hollow sphere with din (inner diameter), dout (outer diameter), and Din
%           'impulsed_vin_d_Dex:   two-compartment with vin (intracellular proton fraction), d (intracelular size), and extracelluar diffusion coefficient Dex
%           'impulsed_vin_d_Dex_Din:    two-compartment with vin, d, Dex, Din
%           'impulsed_vin_d_Dex_Din_betaex:    two-compartment with vin, d, Dex, Din, betaex (extracellular diffusion linear dispersion rate)
% NOTES: 
%       1. 1comptHollow can be used to simulate signals but has not been supported for data fitting yet
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
        out = Signal_d_Din(parms,this) ; 
        out = Signal_d_Din_hollow(parms,this) ; 
        out = Signal_vin_d_Dex(parms, this) ; 
        out = Signal_vin_d_Dex_Din(parms, this) ; 
        out = Signal_vin_d_Dex_Din_betaex(parms, this) ; 
    end

    %% other methods
    methods                
        % ------------------------------- Constructor -------------------------------------
        % IMPULSED should take two inputs: a struct (structure model) and an object (DiffusionPulseSequence)
        function this = IMPULSED(varargin)
            if nargin == 0 , return ; end      % return an empty object
            for n = 1:length(varargin)
                if isa(varargin{n},'mig.DiffusionPulseSequence')
                    this.pulse = varargin{n};
                elseif isstruct(varargin{n})
                    inputStructure = varargin{n};
                else
                    error('Inputs for IMPULSED should be StructureModel, DiffusionPulseSequence')
                end
            end    
            
            %------------------------ setting signal function handle and default fitopts --------------------------
            % !!! Note that the DefaultFitopts.x0, lb, ub should all be column vectors !!!
            switch inputStructure.modelName
                case  '1compt'                                                                                          % [din, Din], 1compartment
                    this.parmsLabels = ["d","Din"] ; this.NumParms = 2 ; 
                    this.FcnSignal = @mig.IMPULSED.Signal_d_Din ;  
                    this.defaultFitopts.x0 = [10, 1.5]' ; this.defaultFitopts.lb = [0 0]' ; this.defaultFitopts.ub = [30 4]' ; 
                    this.defaultStructure.modelName = '1compt' ; this.defaultStructure.geometry = 'sphere' ; 
                case  '1comptHollow'                                                                                          % [din, Din], 1compartment
                    this.parmsLabels = ["din","dout","Din"] ; this.NumParms = 3 ; 
                    this.FcnSignal = @mig.IMPULSED.Signal_d_Din_hollow ;  
                    this.defaultFitopts.x0 = [1,10, 1.5]' ; this.defaultFitopts.lb = [0 0 0]' ; this.defaultFitopts.ub = [5 30 4]' ; 
                    this.defaultStructure.modelName = '1compt' ; this.defaultStructure.geometry = 'hollowSphere' ; 
                case 'impulsed_vin_d_Dex'                                                                                      % [vin, d, Dex]
                    this.parmsLabels = ["vin", "d", "Dex"] ; this.NumParms = 3 ; 
                    this.FcnSignal = @mig.IMPULSED.Signal_vin_d_Dex; 
                    this.defaultFitopts.x0 = [0.5, 10, 1.5]' ;  this.defaultFitopts.lb = [0, 0, 0]' ; this.defaultFitopts.ub = [1, 30 4]' ; 
                    this.defaultStructure.modelName = 'impulsed_vin_d_Dex'  ; this.defaultStructure.Din = 1.56 ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere';
                case 'impulsed_vin_d_Dex_Din'                                                                               % [vin, d, Dex, Din]
                    this.parmsLabels = ["vin", "d", "Dex", "Din"] ; this.NumParms = 4 ; 
                    this.FcnSignal = @mig.IMPULSED.Signal_vin_d_Dex_Din; 
                    this.defaultFitopts.x0 = [0.5, 10, 1.5, 1.5]' ; this.defaultFitopts.lb = [0, 0, 0, 0]' ; this.defaultFitopts.ub = [1, 30, 4,4]' ; 
                    this.defaultStructure.modelName = 'impulsed_vin_d_Dex_Din'  ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere' ; 
                case 'impulsed_vin_d_Dex_Din_betaex'                                                                  % [vin, d, Dex, Din, betaex]       
                    this.parmsLabels = ["vin", "d", "Dex", "Din", "betaex"] ; this.NumParms = 5 ; 
                    this.FcnSignal = @mig.IMPULSED.Signal_vin_d_Dex_Din_betaex ;  
                    this.defaultFitopts.x0 = [0.5, 10, 1.5, 1.5, 3]' ; this.defaultFitopts.lb = [0, 0, 0, 0, 0]' ; this.defaultFitopts.ub = [1, 30, 4, 4, 10]' ; 
                    this.defaultStructure.modelName = 'impulsed_vin_d_Dex_Din_betaex'  ; this.defaultStructure.geometry = 'sphere';
                otherwise
                    error('%s: No matached structure.modelName',mfilename) ; 
            end
            
            % ------------------ update structure -----------------------
            this.structure = this.defaultStructure ; 
            names = fieldnames(inputStructure) ; 
            for n=1:length(names)
                if ~isfield(this.structure, names{n}) , error('%s: The input structure has a field that is not a field of mig.IMPULSED.structure',mfilename) ; end
                this.structure.(names{n}) = inputStructure.(names{n}) ; 
            end

            % -------------------- set other defaultFitopts fields --------------------------------
            this.defaultFitopts.fittingMethod = 'normal' ;          % fittingMethod: {'normal', 'general', 'special'}
            this.defaultFitopts.solverName = 'lsqnonlin' ;          % solver name: {'lsqcurvefit','lsqnonlin','fmincon'}
            this.defaultFitopts.noiseModel = 'none' ;       % {'none','simple','rician','meanrician'}
            this.defaultFitopts.options =optimset('Display','off') ; 
            this.defaultFitopts.ADCfit_bmin = 0.19 ;                % minimum b for ADC fitting
            this.defaultFitopts.ADCfit_bmax = 1.05 ;             % maximum b for ADC fitting
            this.defaultFitopts.NumStarts = 50 ; 
            % fmincon 
            this.defaultFitopts.Aeq = [] ; 
            this.defaultFitopts.beq = [] ; 
            this.defaultFitopts.Aineq = [] ; 
            this.defaultFitopts.bineq = [] ; 
            this.defaultFitopts.nonlcon = [] ; 
            % flags for many options
            this.defaultFitopts.flag.denoise = 'n' ; 
            this.defaultFitopts.flag.degibbs = 'n' ; 
            this.defaultFitopts.flag.deivim = 'y' ; 
            this.defaultFitopts.flag.multistart = 'n' ; 
            this.defaultFitopts.flag.parfor = 'n' ; 
            
            if ~iscolumn(this.defaultFitopts.x0) && ~iscolumn(this.defaultFitopts.lb) && ~iscolumn(this.defaultFitopts.ub)
                error('%s: the DefaultFitopts.x0, lb, ub should all be column vectors',mfilename) ; 
            end
            
        end % End of Constructor
        
    end % End of Methods
    
end % End of class
