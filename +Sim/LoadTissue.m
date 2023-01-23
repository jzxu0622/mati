function structure = LoadTissue(structure, file, flag)
%% load tissue function
% x axis is vertical axis, y axis is horizontal axis


%% get N
structure.N = zeros(1, 3) ; 
str = textscan(file.name, flag.fileFormat, 'Delimiter','_ .') ; 
structure.N(1) = str{2} ; structure.N(2) = str{3} ; structure.N(3) = str{4} ; 
% [str,structure.N(1),structure.N(2),structure.N(3),tt]=strread(filename,filename_format,'delimiter','_.') ; 

%% get dimension of the system
structure.Ndim = sum(structure.N ~= 1) ; 

%% load tissue
switch structure.Ndim
    case 1
        structure.tissue = load(fullfile(file.folder, file.name)) ; 
        % check tissue dimensions 
        if length(structure.tissue) ~= structure.N(1)
            error('tissue size wrong!') ; 
        end
    case 2
        % load file
        structure.tissue = load(fullfile(file.folder, file.name)) ; 
        % check tissue dimensions
        if sum(size(structure.tissue) ~= structure.N(1:2)) > eps
            error('tissue size wrong!') ; 
        end
    case 3
        % load file
        temp2d = load(fullfile(file.folder, file.name)); 
        % convert 2d matrix to 3d
        structure.tissue=zeros(structure.N(1),structure.N(2),structure.N(3)) ; 
        for m=1:structure.N(3)
            structure.tissue(:,:,m) = temp2d(((m-1)*prod(structure.N(1))+1):(m*prod(structure.N(1))),:) ; 
        end
        clear temp2d ; 
        % check tissue dimension
        if sum(size(structure.tissue) ~= structure.N) > eps
            error('tissue size wrong!') ; 
        end
end

    %% convert C-version tissue to Matlab-version tissue
    structure.tissue = structure.tissue + 1 ;      

    %% get # of diffusion compartments
    structure.Ncompt = max(structure.tissue(:)) ; 
    
    %% update D/c/Pm/T2 for all tissue compartments
    structure.c = structure.c(1:structure.Ntype) ; 
    structure.D = structure.D(1:structure.Ntype) ; 
    if structure.flag.useT2 == 'y'
        structure.T2 = structure.T2(1:structure.Ntype) ; 
    end

    % update if Ncompt > Ntype, i.e., same type (cells) but different compartments
    if structure.Ncompt>structure.Ntype
        structure.c = [structure.c , repmat(structure.c(structure.Ntype),[1, structure.Ncompt-structure.Ntype])] ; 
        structure.D = [structure.D , repmat(structure.D(structure.Ntype),[1, structure.Ncompt-structure.Ntype])] ; 
        if structure.flag.useT2 == 'y'
            structure.T2 = [structure.T2 , repmat(structure.T2(structure.Ntype),[1, structure.Ncompt-structure.Ntype])] ; 
        end
    end        

    end
