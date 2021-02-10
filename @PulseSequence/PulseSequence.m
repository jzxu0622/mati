classdef PulseSequence
    %% class definition of PulseSequence, the superclass of all pulse sequence classes
    % NOTE:
    %           1. All inputs should be either a struct or (Nacq, [options]) where options are pairs of names and values. The values should be either a scalar or a vector
    %           2. The output object has all properties as column vectors for calculation convenience
    %
    % ------------------------------------------------------------------------------------------------------------------
    % Author: Junzhong Xu, April 25, 2020
    %
    % ------------------------------------------------------------------------------------------------------------------
    
    
    %% properties
    properties (Constant)
        gamma = 26.75 ;
    end
    
    properties (GetAccess='public', SetAccess='private')
        Nacq ;         % total number of acquisition points
        
        TR {mustBeNumeric, mustBePositive}          % repetition time [ms]
        TE {mustBeNumeric, mustBeNonnegative}          % echo time [ms], could be 0 or positive
        B0 {mustBeNumeric, mustBeNonnegative}          % B0 field strength [Hz]
        df {mustBeNumeric}           % frequency offset [Hz], could be positve or negative
    end
    
   
    %% static methods
    methods (Static)
        function out = cat(varargin)
            % Override cat function by concacenating each field along column
            for narg = 1:nargin
                if ~isa(varargin{narg},'mig.PulseSequence')
                    error('Pulse concatenate takes PulseSequence objects as arguments') ; 
                elseif narg == 1
                    mc = metaclass(varargin{1});
                    is_Constant = [mc.PropertyList.Constant];
                    is_Dependent = [mc.PropertyList.Dependent];
                    all_fields = {mc.PropertyList.Name};
                    fields = all_fields(~is_Constant & ~is_Dependent);
                    for nfield = 1:length(fields)
                        tmp.(fields{nfield}) = varargin{narg}.(fields{nfield});
                    end
                elseif narg > 1
                    for nfield = 1:length(fields)
                        if ~isprop(varargin{narg},fields{nfield})
                            error('Concatenated Pulses must have same property fields');
                        end
                        if strcmp(fields{nfield}, 'Nacq')
                            tmp.Nacq = tmp.Nacq + varargin{narg}.Nacq ; 
                        else
                            if isempty(tmp.(fields{nfield})) ~= isempty(varargin{narg}.(fields{nfield})) % two concatenated fields should either empty or nonempty
                                error('%s: Two concatenated fields %s should be both either empty or nonempty',mfilename, fields{nfield}) ; 
                            end
                            tmp.(fields{nfield}) = [tmp.(fields{nfield}); varargin{narg}.(fields{nfield})];
                        end
                    end
                end
            end
            % ----------------------- output a new object by calling a constractor ---------------------------
            % Note: for convenience, Please add all new subclass constractors here
            if isa(varargin{narg},'mig.DiffusionPulseSequence')
                out = mig.DiffusionPulseSequence(tmp) ; 
            elseif isa(varargin{narg},'mig.SIRqMTPulseSequence') 
            else        % basic PulseSequence
                out = mig.PulseSequence(tmp) ; 
            end
        end
    
        function out = disp(this)
            mc = metaclass(this) ; 
            fields = {mc.PropertyList.Name} ; 
            for n=1:length(fields)
                out.(fields{n}) = this.(fields{n})' ;  
            end
            mig.Tools.cprintf('_blue','Note: PulseSequence properties should be column vectors. They are shown as row vectors here for visualization purpose only\n') ; 
        end
    
    end     % end of static methods
    
    %% methods
    methods
        function this = set.Nacq(this,Nacq) ,             this.Nacq = Nacq ;         end
        function this = set.TE(this,input) ,            this.TE = assignValue(this,input,1) ;         end
        function this = set.TR(this,input) ,             this.TR = assignValue(this,input,1) ;         end
        function this = set.B0(this,input) ,             this.B0 = assignValue(this,input,1) ;         end
        function this = set.df(this,input) ,             this.df = assignValue(this,input,1) ;         end

        % %%%%%%%%% constructor %%%%%%%%%%%%%%%%%
        function this = PulseSequence(varargin)
            % -------------- read all inputs ------------------
            switch nargin
                case 0
                    return ; 
                case 1      % 1 input should be a struct
                    if ~isstruct(varargin{1}) , error('%s: If there is one input, it should be a struct',mfilename) ; end
                    names = fieldnames(varargin{1});
                    values = cell(size(names));
                    if ~ismember('Nacq',names) , error('%s: Nacq is a required field in the input struct',mfilename) ; end
                    for i = 1:length(names)
                        values{i} = varargin{1}.(names{i});
                        if strcmp(names{i},'Nacq'), Nacq = values{i}; end
                    end
                case 2
                    error('%s: The input should be a struct or (Nacq, ''property'', value, ''property'', value, ...])',mfilename) ; 
                otherwise        % nargin >= 3, construction using pairs of ['property',value]
                    Nacq = varargin{1};     % Nacq always the first input arg
                    names = varargin(2:2:end);
                    values = varargin(3:2:end);
                    if length(unique(names)) ~= length(names), error('%s: The input names have duplicates',mfilename) ; end
                    if length(names) ~= length(values) , error('%s: The input should have pairs of ''property'' and value])',mfilename) ; end
            end
            % --------------- construct object ----------------
            if (rem(Nacq,1) ~= 0) || (Nacq < 0), error('%s: Nacq should be a positive integer', mfilename) ; end
            this.Nacq = Nacq ; 
            % assign independent properties
            for narg=1:length(names)
                if ~isstring(names{narg}) && ~ischar(names{narg})
                    error('%s: The input should be a struct or has a form of (Nacq,[options]) where [options] is list of name, value pairs', mfilename)
                end
                if isprop(this,names{narg}) % If name is a property
                    if ~strcmp(names{narg},'gamma') && ~strcmp(names{narg},'Nacq') % Can't set constant property
                        this.(names{narg}) = assignValue(this, values{narg},1);
                    end
                else
                    error('%s: Unknown input keyword %s',mfilename, names{narg});
                end
            end
            % --------------------------------- check inputs -------------------------------------
            % ----------------------- calculate dependent properties -----------------------------
        
        end
        
        %%%%%%%%%%%%%% function to assign values %%%%%%%%%%%%%%%
%         function output = assignValue(this, input)
%         % set parameter as values. input can be a scalar or a vector. output should be a column vector
%             if ischar(input) , input = convertCharsToStrings(input) ; end
%             if isscalar(input)
%                 output = repmat(input, [this.Nacq,1]) ;
%             elseif isempty(input)
%                 output = [];
%             else
%                 if ~isvector(input), error('%s: The input value should be a scalar or a vector',mfilename); end
%                 if length(input) ~= this.Nacq, error('%s: All input pulse parameters should be a scalar or a vector with a length of Nacq=%d',mfilename,this.Nacq); end
%                 if iscolumn(input)  % make sure input is a column vector
%                     output=input ;
%                 else
%                     output = input' ;
%                 end
%             end
%         end
        
        function output = assignValue(this, input, Ncol)
            if ischar(input) , input = convertCharsToStrings(input) ; end
            if isempty(input)
                output = [] ; 
            else
                switch Ncol 
                    case 1         % set parameters of non-directional
                    if isscalar(input)
                        output = repmat(input, [this.Nacq,1]) ;
                    else
                        if ~isvector(input), error('%s: The input value should be a scalar or a vector',mfilename); end
                        if length(input) ~= this.Nacq, error('%s: All input pulse parameters should be a scalar or a vector with a length of Nacq=%d',mfilename,this.Nacq); end
                        if iscolumn(input)  % make sure input is a column vector
                            output=input ;
                        else
                            output = input' ;
                        end
                    end
                    case 3        % directional parameters
                        if isvector(input)
                            if iscolumn(input), input = input' ; end
                            output = repmat(input,[this.Nacq,1]) ; 
                        else    %matrix
                            output = input ; 
                        end
                    otherwise
                        error('%s: The input Ncol should be either 1 for all parameters or 3 for e.g., directions',mfilename) ; 
                end
            end
        end % end of assignValue
            
    end % end of methods
            
end % end of class PulseSequence
