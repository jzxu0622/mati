classdef constADC < mati.SignalModel
%% class definition of constADC  
% the same ADC for all DWI measurements 
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
        out = Signal_constADC_mono(parms, this) ;
        out = Signal_constADC_biexp(parms, this) ;        
    end

    %% other methods
    methods                
        % ------------------------------- Constructor -------------------------------------
        % constADC should take two inputs: a struct (structure model) and an object (DiffusionPulseSequence)
        function this = constADC(varargin)
            if nargin == 0 , return ; end      % return an empty object
            for n = 1:length(varargin)
                if isa(varargin{n},'mati.DiffusionPulseSequence')
                    this.pulse = varargin{n};
                elseif isstruct(varargin{n})
                    inputStructure = varargin{n};
                else
                    error('Inputs for constADC should be StructureModel, DiffusionPulseSequence')
                end
            end    
            
            %------------------------ setting signal function handle and default fitopts --------------------------
            % !!! Note that the DefaultFitopts.x0, lb, ub should all be column vectors !!!
            switch inputStructure.modelName
                case 'constADC_mono'
                    this.parmsLabels = ["S0", "ADC"] ; 
                    this.FcnSignal = @mati.constADC.Signal_constADC_mono; 
                    this.defaultFitopts.x0 = [1, 1]' ; this.defaultFitopts.lb = [0, 0]' ; this.defaultFitopts.ub = [1, 3]' ; 
                    this.defaultStructure.modelName = 'constADC_mono'  ; 
                case 'constADC_biexp'
                    this.parmsLabels = ["S0", "f", "ADC1", "ADC2"] ; 
                    this.FcnSignal = @mati.constADC.Signal_constADC_biexp; 
                    this.defaultFitopts.x0 = [1, 0.5, 1, 1]' ; this.defaultFitopts.lb = [0, 0, 0, 0]' ; this.defaultFitopts.ub = [1, 1, 3, 3]' ; 
                    this.defaultStructure.modelName = 'constADC_biexp'  ;  
                otherwise
                    error('%s: No matached structure.modelName',mfilename) ; 
            end
            this.NumParms = length(this.parmsLabels);            
            % ------------------ update structure -----------------------
            this.structure = this.defaultStructure ; 
            names = fieldnames(inputStructure) ; 
            for n=1:length(names)
                if ~isfield(this.structure, names{n}) , error('%s: The input structure has a field that is not a field of mati.constADC.structure',mfilename) ; end
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
