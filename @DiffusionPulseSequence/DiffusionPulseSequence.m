classdef DiffusionPulseSequence < mig.PulseSequence
%% class definition of DiffusionPulseSequence
%
% -------------------------------------------------------------------------------------------------------------------------
% NOTES:
%           1. All inputs should be in pairs ('property', value), and the values should be either a scalar or a vector
%           2. the output object has all properties as column vectors
% -------------------------------------------------------------------------------------------------------------------------
% Author: Junzhong Xu, April 25, 2020
%
% -------------------------------------------------------------------------------------------------------------------------
    
    
    %% properties
    properties (GetAccess='public', SetAccess='private')
        delta {mustBeNumeric, mustBePositive}       % gradient duration
        Delta {mustBeNumeric, mustBePositive}       % gradient separation
        shape {mustBeMember(shape,{'pgse','tpgse','cos','tcos','sin','tsin'})} = 'tpgse' % gradient shape, a single string or a string array. The default is 'pgse'
        b {mustBeNumeric, mustBeNonnegative}       % b value
        G {mustBeNumeric, mustBeNonnegative}      % maximum gradient strength
        n {mustBeNumeric,mustBeInteger}                 % number of gradient oscillating cycles
        trise {mustBeNumeric, mustBeNonnegative}  % gradient rise time, default = 0
        gdir {mustBeNumeric}             % [rx,ry,rz] gradient direction
    end
    properties (Dependent)        % dependent properties
        f ;                % gradient frequency
        w ;               % angular gradient frequency
        T ;               % gradient shape period
        tp ;              % gradient plateau time
        tdiff ;           % effective diffusion time
    end
    
   
    %% methods
    methods
        function this = set.delta(this,input) ,             this.delta = assignValue(this,input,1) ;         end
        function this = set.Delta(this,input) ,             this.Delta = assignValue(this,input,1) ;         end
        function this = set.b(this,input) ,             this.b = assignValue(this,input,1) ;         end
        function this = set.G(this,input) ,             this.G = assignValue(this,input,1) ;         end
        function this = set.n(this,input) ,             this.n = assignValue(this,input,1) ;         end
        function this = set.trise(this,input) ,             this.trise = assignValue(this,input,1) ;         end
        function this = set.shape(this,input) ,             this.shape = assignValue(this,input,1) ;         end
        % Set gdir
        function this = set.gdir(this,input) ,            this.gdir = assignValue(this, input,3) ; end
        
        % ----------------- Dependent Properties ------------------------
        function out = get.f(this)
            out = (this.n > 0).*this.n ./ this.delta ...
                + (this.n == 0) .* 0 ; % (this.n == 0).*(1/4./this.tdiff);
        end
        function out = get.w(this) ,             out = 2*pi*this.f ;         end
        function out = get.T(this) ,             out = 1./this.f ;             out(this.n==0) = 0 ;         end
        % get tp
        function out = get.tp(this)
            if ~ismember("pgse",this.shape) && ~ismember("tpgse",this.shape) && ~ismember("tcos",this.shape) % only pgse, tpgse, and tcos require tp
                out.tp = [] ;             
            elseif isempty(this.trise)
                error('%s: The shape has either pgse, tpgse, or tcos so trise must be set',mfilename) ; 
            else
                out = zeros([this.Nacq,1]) ; 
                ind = (this.shape=="pgse") ;   if ~isempty(ind) ,  out(ind) = this.delta(ind) ; end
                ind = (this.shape=="tpgse") ;   if ~isempty(ind) ,  out(ind) = this.delta(ind) - 2*this.trise(ind) ; end
                ind = (this.shape=="tcos") ;   if ~isempty(ind) ,  out(ind) =  (this.delta(ind) - (6*this.n(ind)+1).*this.trise(ind))/4./this.n(ind); end
            end
            if any(out<0), error('%s: tp cannot be <0',mfilename) ; end
        end
        % get tdiff
        function out = get.tdiff(this)
            out = this.Delta - this.delta/3 ;       % default is pgse
            ind = (this.shape=="sin") ;   if ~isempty(ind) ,  out(ind) = 3/4./(this.n(ind) ./ this.delta(ind)) ; end
            ind = (this.shape=="cos") ;   if ~isempty(ind) ,  out(ind) = 1/4./(this.n(ind) ./ this.delta(ind)) ; end
            ind = (this.shape=="tcos") ;   if ~isempty(ind) ,  out(ind) = 1/4./(this.n(ind) ./ this.delta(ind)) ; end
        end
        
        % %%%%%%%%% Constructor %%%%%%%%%%%%%%%%%
        function this = DiffusionPulseSequence(varargin)
            % -------------- read all inputs ------------------
            superStruct.Nacq = 0 ; subStruct = [] ; 
            switch nargin
                case 0      % empty object
                    names = [] ; 
                case 1      % One struct input
                    if ~isstruct(varargin{1}) , error('%s: If there is one input, it should be a struct',mfilename) ; end
                    names = fieldnames(varargin{1});
                    values = cell(size(names));
                    for i = 1:length(names)
                        values{i} = varargin{1}.(names{i});
                        if strcmp(names{i},'Nacq'), superStruct.Nacq = values{i}; end
                    end
                case 2      % pairs of ['propery',value]
                    error('%s: The input should be a struct or (Nacq, ''property'', value, ''property'', value, ...])',mfilename) ; 
                otherwise        % nargin >= 3, construction using pairs of ['property',value]
                    superStruct.Nacq = varargin{1};     % Nacq always the first input arg
                    names = varargin(2:2:end);
                    values = varargin(3:2:end);
                    if size(names) ~= size(values) , error('%s: The input should have pairs of [''property'',value])',mfilename) ; end
                    for narg=1:length(names)
                        if ~isstring(names{narg}) && ~ischar(names{narg}) , error('%s: The input should be (struct) or (Nacq,[""property"",value,...]', mfilename) ; end
                    end
            end
            % --------------- construct object ----------------------
            superPropertyNames = properties(mig.PulseSequence) ; % find properties of super(parent)-class that need to call superclass constructor
            for n=1:length(names)
                if ismember(names{n},superPropertyNames)
                    superStruct.(names{n}) = values{n} ;  
                else
                    subStruct.(names{n}) = values{n} ; 
                end
            end
            % construct a superclass
            this = this@mig.PulseSequence(superStruct) ; 
            % continue to construct a subclass
            subPropertyNames = properties(this) ; 
            if isempty(subStruct)  % no input so return an empty object
                return
            end
            % set all diffusion properties and check values
            names = fieldnames(subStruct) ; 
            for n=1:length(names)
                if ~ismember(names{n}, subPropertyNames), error('%s: The input property name is not recognized',mfilename) ; end
                this.(names{n}) = subStruct.(names{n}) ; 
            end
            
            % ------------------------ check inputs ---------------------------------
            if isempty(this.delta) , error('%s: The input gradient duration delta is not set',mfilename) ; end
            if isempty(this.Delta) 
                if ismember('pgse',this.shape) || ismember('tpgse',this.shape)
                    error('%s: The input gradient separation Delta should be set for pgse or tpgse sequences',mfilename) ; 
                end
            end
            % Check if n matches shape
            ind_pgse =( this.shape == 'pgse' | this.shape == 'tpgse') ;
            if isempty(this.n)
                this.n = zeros(size(this.shape)) ; 
                if sum(ind_pgse) ~= length(this.shape)
                    error('%s: The number of oscillations ''n'' must be defined for OGSE sequences',mfilename) ; 
                end
            elseif any(this.n(ind_pgse) > 0)
                error('%s: n should be zero for PGSE', mfilename)
            end
            % check trise
            if any(this.shape=='tpgse' | this.shape=='tcos') && isempty(this.trise), error('%s: The trise must be set for tpgse or tcos shapes',mfilename) ; end 
            % check b and G
            if isempty(this.b) && isempty(this.G), error('%s: either b or G should be set',mfilename); end
            if ~isempty(this.b)
                if ~isempty(this.G), warning('%s: Both b and G were set. The input G will be ignored !',mfilename); end
                flag.b_fixed = 'y' ;
                this.G = zeros(size(this.b)) ;
            else
                flag.b_fixed = 'n' ;
                this.b = zeros(size(this.G)) ;
            end
            % check g
            if ~isempty(this.gdir)
                if size(this.gdir,1) ~= this.Nacq || size(this.gdir,2) ~= 3         % not a Nacq x 3 matrix
                    if size(this.gdir,1) == 3 || size(this.gdir,2) == this.Nacq     % a 3 x Nacq matrix
                        warning('%s: gdir appears a 3 x Nacq matrix. change it to a Nacq x 3 matrix',mfilename) ; 
                        this.gdir = this.gdir' ; 
                    else
                        error('%s: The input gdir should be a vector or a Nacq x 3 matrix',mfilename) ; 
                    end
                end
                % normalize direction vectors
                this.gdir = this.gdir ./ repmat(sqrt(sum(this.gdir.^2,2)),[1,3]) ; 
            end
            % ------------------------------ calculate dependent properties -----------------------------
            ind = (lower(this.shape)=='tpgse') ;     % trapzezoidal PGSE
            if any(ind)
                if flag.b_fixed == 'y'    % calculate G using the input b
                    this.G(ind) = sqrt(this.b(ind)./((this.trise(ind)+this.tp(ind)).^2.*(this.Delta(ind)-(this.trise(ind)+this.tp(ind))/3)+this.trise(ind).^3/30 - (this.trise(ind)+this.tp(ind)).*this.trise(ind).^2/6)) / this.gamma ;
                else % calculate b using the input G. Note that the previous check makes sure b and G cannot be empty simultaneously
                    this.b(ind) = (this.gamma * this.G(ind)).^2 .* ((this.trise(ind)+this.tp(ind)).^2 .*(this.Delta(ind)-(this.trise(ind)+this.tp(ind))/3)+this.trise(ind).^3/30 - (this.trise(ind)+this.tp(ind)).*this.trise(ind).^2/6) ;
                end
            end
            ind = (lower(this.shape)=='tcos') ;       % trapezoidal cos
            if any(ind)
                if flag.b_fixed == 'y'
                    this.G(ind) = sqrt(this.b(ind)./((91*this.n(ind).*this.trise(ind).^3)/15 + (8*this.n(ind).*this.tp(ind).^3)/3 + this.trise(ind).^3/30 + 12*this.n(ind).*this.trise(ind).*this.tp(ind).^2 + (46*this.n(ind).*this.trise(ind).^2.*this.tp(ind))/3)) / this.gamma ;
                else
                    this.b(ind) = (this.gamma*this.G(ind)).^2.*((91*this.n(ind).*this.trise(ind).^3)/15 + (8*this.n(ind).*this.tp(ind).^3)/3 + this.trise(ind).^3/30 + 12*this.n(ind).*this.trise(ind).*this.tp(ind).^2 + (46*this.n(ind).*this.trise(ind).^2.*this.tp(ind))/3) ;
                end
            end
            ind = (lower(this.shape)=='pgse') ;      % rectangular PGSE
            if any(ind)
                if flag.b_fixed == 'y'
                    this.G(ind) = sqrt(this.b(ind)./(this.Delta(ind)-this.delta(ind)/3)) / this.gamma ./ this.delta(ind) ;
                else
                    this.b(ind) = (this.gamma .* this.G(ind) .* this.delta(ind)).^2 .* (this.Delta(ind) - this.delta(ind)/3) ;
                end
            end
            ind = (lower(this.shape)=='cos') ;        % cos
            if any(ind)
                if flag.b_fixed == 'y'
                    this.G(ind) = sqrt(4*this.b(ind)./this.delta(ind).^3).*this.n(ind)*pi/this.gamma ;
                else
                    this.b(ind) = 1./4 * (this.gamma*this.G(ind)/pi./this.n(ind)).^2 .* this.delta(ind).^3 ;
                end
            end
            ind = (lower(this.shape)=='sin') ;         % sin
            if any(ind)
                if flag.b_fixed == 'y'
                    this.G(ind) = sqrt(this.b(ind)*4/3./this.delta(ind).^3) * pi .* this.n(ind) / this.gamma ;
                else
                    this.b(ind) = 3/4*(this.gamma*this.G(ind)/pi./ this.n(ind)).^2 .* this.delta(ind).^3 ;
                end
            end
            ind = (lower(this.shape)=='tsin') ;          % trapezoidal sin
            if any(ind) , error('%s: The tsin shape is not supported yet',mfilename) ; end

        end % end of constructor
        
    end % end of methods
    
end % end of class DiffusionPulseSequence
