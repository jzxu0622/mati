classdef FitPars
%% class defination of FitPars to save fitting parameters
%
% -------------------------------------------------------------------------------------------------------------------------

    properties
        model
        fitopts
    end
    methods (Static)
        out = ObjectiveFcn(this, parms, ydata, sigma) ; 
    end
    
    
    methods        
        function this = set.fitopts(this,fitopts)
            this.fitopts = fitopts;
        end
        function this = set.model(this,model)
            this.model = model;
        end
        
        function this = FitPars(varargin)
            switch nargin
                case 0
                    return 
                case 1
                    if ~isa(varargin{1}, 'mig.SignalModel'), error('%s: FitPars should take a SignalModel object', mfilename) ; end
                    this.model = varargin{1} ; 
                    this.fitopts = this.model.defaultFitopts ; 
                case 2
                    for n = 1:nargin        % fitopts structure and SignalModel object can be in any order
                        if isa(varargin{n},'mig.SignalModel')
                            this.model = varargin{n} ; 
                        elseif isstruct(varargin{n})
                            inputFitopts = varargin{n} ;              % take input fitopts structure
                        else
                            error('FitPars can only take a fitopts structure and SignalModel object inputs')
                        end
                    end
                    % use input fitopts to override defaultFitopts fields
                    this.fitopts = this.model.defaultFitopts ;  % place holder for all default fitopts
                    default_names = fieldnames(this.model.defaultFitopts) ; 
                    inputNames = fieldnames(inputFitopts) ; 
                    for n = 1:length(inputNames)
                        if strcmp(inputNames{n}, 'flag')
                            defaultFlagNames = fieldnames(this.model.defaultFitopts.flag) ; 
                            inputFlagNames = fieldnames(inputFitopts.flag) ; 
                            for m = 1:length(inputFlagNames)
                                if ~ismember(inputFlagNames{m}, defaultFlagNames) 
                                    warning on ; warning('%s: The input fitopts.flag.%s is not a new field, not in defaultFitopts.flag',mfilename, inputFlagNames{m}); warning off ; end
                                this.fitopts.flag.(inputFlagNames{m}) = inputFitopts.flag.(inputFlagNames{m}) ; 
                            end
                        else
                            if ~ismember(inputNames{n}, default_names), warning on ; warning('%s: The input fitopts.%s is a new field that is not in mig.IMPULSED.defaultFitopts',mfilename, inputNames{n}); warning off ; end
                            this.fitopts.(inputNames{n}) = inputFitopts.(inputNames{n}) ; 
                        end
                    end
    
                otherwise
                    error('%s: Input should be either a SignalModel object,either with or without a fitopts structure', mfilename) ; 
            end % End of switch nargin
            
        end % End of constructor
        
    end % End of Methods
    
end % End of class
