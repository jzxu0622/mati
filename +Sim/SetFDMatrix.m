function structure = SetFDMatrix(structure) 
%% set transfer matrix
% get the transverse matrices for both 2D and 3D
% output: A -- constant elements of B
%         Ax1 -- upper(negative) B elements dependent on the x gradient, periodic BC
%         Ax2 -- lower(positive) B elements dependent on the x gradient, periodic BC
%         Ay1 -- upper(negative) B elements dependent on the y gradient, periodic BC
%         Ay2 -- lower(positive) B elements dependent on the y gradient, periodic BC
%         Az1 -- upper(negative) B elements dependent on the z gradient, periodic BC
%         Az2 -- lower(positive) B elements dependent on the z gradient, periodic BC
% subfunction: jp -- determine the jump probability 
% *************************************************************************
% row: x
% column: y
% layer: z
% *************************************************************************


    %% initialize matrices
    % total size of tissue matrix
    structure.size = prod(structure.N) ; 
    
    % set C (concentration) and D (diffusivity) matrices for computing convenience
    C = zeros(size(structure.tissue)) ; D = C ; 
    for n = 1:structure.Ncompt
        C(structure.tissue==n) = structure.c(n) ; 
        D(structure.tissue==n) = structure.D(n) ; 
    end

    % index matrix
    IND = reshape(1:structure.size, structure.N) ; 
    % coordination matrix
    COORD = cell([3,1]) ; X = 1:structure.N(1) ; Y = 1:structure.N(2) ; Z = 1:structure.N(3) ; [COORD{1}, COORD{2}, COORD{3}] = ndgrid(X, Y, Z) ; 
    Adiag = ones([structure.size,1]) ; 

    % initialization of FD matrices
    FDMatrix = cell([structure.Ndim,1]) ; FDMatrixEdge = FDMatrix ; 
    for ndim = 1:3  % x,y,z directions
        if structure.Ndim<ndim
            FDMatrix{ndim}{1} = sparse(structure.size, structure.size) ;             FDMatrix{ndim}{2} = FDMatrix{ndim}{1} ; 
            FDMatrixEdge{ndim}{1} = FDMatrix{ndim}{1} ; FDMatrixEdge{ndim}{2} = FDMatrix{ndim}{1} ; 
            continue ; 
        else
            for n =1:2 % diffusion towards negative axis (n=1) or positive axis (n=2)
                if n==1
                    ndir = 1 ; NEdge = structure.N(ndim) ; 
                else
                    ndir = -1 ; NEdge =  1; 
                end

                % circular shift matrices 
                INDShift = circshift(IND, ndir, ndim) ;  DShift = circshift(D, ndir, ndim) ; CShift = circshift(C,ndir,ndim) ; TissueShift = circshift(structure.tissue,ndir,ndim) ; 
                JPC = D(:) * structure.dt / structure.dr(ndim)^2 ; % default jpc inside the same compartment
                if structure.Pm > 100   % completel permeable
                    JPCex = 2 * structure.dt / structure.dr(ndim)^2 *D(:) .*DShift(:) ./ (D(:).*C(:) + DShift(:).*CShift(:)) ; 
                else
                    JPCex = 2 * structure.dt /structure.dr(ndim) .*D(:) .*DShift(:) * structure.Pm *structure.cfree ./ (structure.Pm*structure.dr(ndim)*structure.cfree*D(:).*C(:) + 2*D(:).*C(:).*DShift(:).*CShift(:) + structure.Pm*structure.dr(ndim)*structure.cfree*DShift(:).*CShift(:)) ;
                end
                % correction for jpc=nan from c=0 to c=0 compartments
                JPCex(isnan(JPCex(:))) = 0 ; 
                % update all JPC elements with Pm
                indPm = structure.tissue(:) ~= TissueShift(:) ; 
                JPC(indPm) = JPCex(indPm) ; 

                % set indices of matrix elements
                ind = circshift(COORD{ndim}, ndir,ndim) ; ind = ind(:) ; ind = (ind~=NEdge) ; indEdge = ~ind ; 
                
                % construct FD matrices using graphy
                G = digraph(IND(:), INDShift(:)) ; 
                G.Edges.Weight = CShift(:) .* JPC .* ind ; 
                FDMatrix{ndim}{n} = adjacency(G, 'weighted') ; 
                Adiag = Adiag - C(:) .* JPC .* ind ; 

                % update edge elements w/ periodic boundary condition
                if strcmpi(structure.BC, 'rpbc') 
                    G = digraph(IND(:), INDShift) ; 
                    G.Edges.Weight = CShift(:) .* JPC .* indEdge ; 
                    FDMatrixEdge{ndim}{n} = adjacency(G, 'weighted') ; 
                    Adiag = Adiag - C(:) .* JPC .* indEdge ; 
                end
            end
        end
    end
    
    %% finalize FD matrices
    structure.FDMatrix.A = sparse(1:structure.size, 1:structure.size, Adiag) + FDMatrix{1}{1}  + FDMatrix{1}{2} + FDMatrix{2}{1}  + FDMatrix{2}{2} + FDMatrix{3}{1}  + FDMatrix{3}{2} ; 
    if strcmpi(structure.BC, 'rpbc')
        structure.FDMatrix.Ax1 = FDMatrixEdge{1}{1} ; 
        structure.FDMatrix.Ax2 = FDMatrixEdge{1}{2} ; 
        structure.FDMatrix.Ay1 = FDMatrixEdge{2}{1} ; 
        structure.FDMatrix.Ay2 = FDMatrixEdge{2}{2} ; 
        structure.FDMatrix.Az1 = FDMatrixEdge{3}{1} ; 
        structure.FDMatrix.Az2 = FDMatrixEdge{3}{2} ; 
    else    % bounded boundary condition. Set all edge transit matrices to empty sparse matrices
        structure.FDMatrix.Ax1 = sparse(structure.size, structure.size)  ; 
        structure.FDMatrix.Ax2 = sparse(structure.size, structure.size)  ; 
        structure.FDMatrix.Ay1 = sparse(structure.size, structure.size)  ; 
        structure.FDMatrix.Ay2 = sparse(structure.size, structure.size)  ; 
        structure.FDMatrix.Az1 = sparse(structure.size, structure.size)  ; 
        structure.FDMatrix.Az2 = sparse(structure.size, structure.size)  ;         
    end

end


