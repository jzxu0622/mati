classdef JOINT < mati.SignalModel
%% class definition of JOINT  
% JOINT can take different structure models. For OGSE signals, water exchange effects are ignored. Fore PGSE signals, water exchange effects are included using Karger model.
% Model choices: 
%           'joint_vin_d_Dex_Din_kin:   two-compartment with vin (intracellular proton fraction), d (intracelular size), Dex (extracelluar diffusion coefficient), and kin (water exchange rate)
%           'joint_vin_d_Dex_kin:    two-compartment with vin, d, kin and a fixed Din      
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
        out = Signal_joint_vin_d_Dex_Din_kin(parms, this) ;
        out = Signal_joint_vin_d_Dex_kin(parms, this) ;
        [ M_in, M_ex ] = generalSolution_blochEQ( Ain,Aex,kin,kex,Min0,Mex0,t ) ; 
    end

    %% other methods
    methods                
        % ------------------------------- Constructor -------------------------------------
        % JOINT should take two inputs: a struct (structure model) and an object (DiffusionPulseSequence)
        function this = JOINT(varargin)
            if nargin == 0 , return ; end      % return an empty object
            for n = 1:length(varargin)
                if isa(varargin{n},'mati.DiffusionPulseSequence')
                    this.pulse = varargin{n};
                elseif isstruct(varargin{n})
                    inputStructure = varargin{n};
                else
                    error('Inputs for JOINT should be StructureModel, DiffusionPulseSequence')
                end
            end    
            
            %------------------------ setting signal function handle and default fitopts --------------------------
            % !!! Note that the DefaultFitopts.x0, lb, ub should all be column vectors !!!
            switch inputStructure.modelName
                case 'joint_vin_d_Dex_Din_kin'
                    this.parmsLabels = ["vin", "d", "Dex", "Din", "kin"] ; 
                    this.FcnSignal = @mati.JOINT.Signal_joint_vin_d_Dex_Din_kin; 
                    this.defaultFitopts.x0 = [0.5, 10, 1.5, 1.56, 1/200]' ; this.defaultFitopts.lb = [0, 0, 0, 0, 1/1000]' ; this.defaultFitopts.ub = [1, 30, 3,3,1/10]' ; 
                    this.defaultStructure.modelName = 'joint_vin_d_Dex_Din_kin'  ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere' ; 
                case 'joint_vin_d_Dex_kin'
                    this.parmsLabels = ["vin", "d", "Dex", "kin"] ; 
                    this.FcnSignal = @mati.JOINT.Signal_joint_vin_d_Dex_kin; 
                    this.defaultFitopts.x0 = [0.5, 10, 1.5, 1/200]' ; this.defaultFitopts.lb = [0, 0, 0, 1/1000]' ; this.defaultFitopts.ub = [1, 30, 3,1/10]' ; 
                    this.defaultStructure.modelName = 'joint_vin_d_Dex_kin'  ; this.defaultStructure.Din = 1.56 ; this.defaultStructure.betaex = 0 ; this.defaultStructure.geometry = 'sphere' ; 
                otherwise
                    error('%s: No matached structure.modelName',mfilename) ; 
            end
            this.NumParms = length(this.parmsLabels);            
            % ------------------ update structure -----------------------
            this.structure = this.defaultStructure ; 
            names = fieldnames(inputStructure) ; 
            for n=1:length(names)
                if ~isfield(this.structure, names{n}) , error('%s: The input structure has a field that is not a field of mati.JOINT.structure',mfilename) ; end
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
