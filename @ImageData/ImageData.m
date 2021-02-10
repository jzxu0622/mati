classdef ImageData
%% class definition of ImageData which stores image data and some parameters
%
% -------------------------------------------------------------------------------------------------------------------------

    properties        
        img {mustBeNumeric}
        sigma {mustBeNumeric}
        mask {mustBeNumericOrLogical}
    end
    
    properties (Dependent)
        Nx
        Ny
        Nz
        Nt
    end
    
    methods
        function out = get.Nx(this)
            if ~isempty(this.img) , out = size(this.img,1) ; else , out = [] ; end
        end
        function out = get.Ny(this)
            if ~isempty(this.img) , out = size(this.img,2) ; else , out = [] ; end
        end
        function out = get.Nz(this)
            if ~isempty(this.img) , out = size(this.img,3) ; else , out = [] ; end
        end
        function out = get.Nt(this)
            if ~isempty(this.img) , out = size(this.img,4) ; else , out = [] ; end
        end
        
        %%%%%%%%%%%%%% Constructor %%%%%%%%%%%%%%%%%%
        function this = ImageData(varargin)
            if nargin == 0
                return
            else        % at least one input of image matrix
                this.img = varargin{1} ; 
                if ~isnumeric(this.img) , error('%s: The input image data must be a Nx x Ny x Nz x Nt numeric matrix',mfilename) ; end
                if ndims(this.img) > 4 , error('%s: The input matrix must be a Nx x Ny x Nz x Nt numeric matrix',mfilename) ; end
                [Nx,Ny,Nz,Nt] = size(this.img) ; 
                this.sigma = [] ; 
                this.mask = ones([Nx,Ny,Nz])>0 ; 
                if nargin > 1 && ~isempty(varargin{2})
                    this.sigma = varargin{2} ; 
                    if ~isnumeric(this.sigma) , error('%s: The input noise sigma must be empty, a scalar, or a Nx x Ny x Nz  numeric matrix',mfilename) ; end
                    if length(size(this.sigma)) > 3 , error('%s: The input noise sigma must be a Nx x Ny x Nz matrix',mfilename) ; end
                    if isscalar(this.sigma) , this.sigma = repmat(this.sigma, [Nx, Ny, Nz]) ; end
                    if size(this.img,1)~=size(this.sigma,1) || size(this.img,2)~=size(this.sigma,2) || size(this.img,3)~=size(this.sigma,3)
                        error('%s: The first three dimensions of image data and noise sigma are not equal',mfilename) ; 
                    end
                end
                if nargin == 3
                    this.mask = varargin{3}>0 ; 
                    if isempty(this.mask) , error('%s: The input mask cannot be empty',mfilename) ; end 
                    if ~isnumeric(this.mask) && ~islogical(this.mask), error('%s: The input mask must be a Nx x Ny x Nz  numeric or logical matrix',mfilename) ; end
                    if length(size(this.mask)) > 3 , error('%s: The input mask must be a Nx x Ny x Nz matrix',mfilename) ; end
                    if size(this.img,1)~=size(this.mask,1) || size(this.img,2)~=size(this.mask,2) || size(this.img,3)~=size(this.mask,3)
                        error('%s: The first three dimensions of image data and mask are not equal',mfilename) ; 
                    end
                end
                if nargin > 3
                    error('%s: Data object currently only takes three inputs: [image data, noise sigma, and mask]',mfilename) ; 
                end
            end
            
        end % End of Constructor
        
    end % End of methods
    
end % End of class definiation

