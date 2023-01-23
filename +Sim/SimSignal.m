function M = SimSignal(structure, pulse, flag)
%% Simulate DWI signal for a certain structure, a pulse sequence 
% OUTPUTS
%       M: a column vector of magnetizations at each simulated point/voxel
%
% ----------------------------------------------------------------------------------------

    %% initialize magnetization as 1D column vector, w/ a length of structure.size
    % \pi/2 rf pulse, spins are flipped onto y-axis
    M = 0 ; 
    for n=1:structure.Ncompt
        M = M + 1i * structure.c(n) .* (structure.tissue(:)==n) ;
    end
       
    
    %% if use GPU, send the vector and FD matrix to GPU
    if flag.useGPU == 'y'
        M = gpuArray(M) ; 

        structure.coordRange = gpuArray(structure.coordRange) ; 
    
        structure.FDMatrix.A = gpuArray(structure.FDMatrix.A) ; 
        structure.FDMatrix.Ax1 = gpuArray(structure.FDMatrix.Ax1) ; 
        structure.FDMatrix.Ax2 = gpuArray(structure.FDMatrix.Ax2) ; 
        structure.FDMatrix.Ay1 = gpuArray(structure.FDMatrix.Ay1) ; 
        structure.FDMatrix.Ay2 = gpuArray(structure.FDMatrix.Ay2) ; 
        structure.FDMatrix.Az1 = gpuArray(structure.FDMatrix.Az1) ; 
        structure.FDMatrix.Az2 = gpuArray(structure.FDMatrix.Az2) ; 
    end
    

    %% Simulation
    for nt=1:length(pulse.grad)
        % diffusion gradient induced phase change
        if abs(pulse.grad(nt)) > eps || nt == 1
            M = M .* exp (- 1i .* pulse.gamma .* pulse.dt .* structure.coordRange .* pulse.grad(nt)) ;  
        end
            
        % updating M in each time step
        if strcmpi(structure.BC,'rpbc')                 % update B if periodic boundary condition
            M = structure.FDMatrix.A * M ... 
                + structure.FDMatrix.Ax1 .* exp( structure.rpbcFactorX .* pulse.gIntegral(nt)) * M ...
                + structure.FDMatrix.Ax2 .* exp(-structure.rpbcFactorX .* pulse.gIntegral(nt)) * M ...
                + structure.FDMatrix.Ay1 .* exp( structure.rpbcFactorY .* pulse.gIntegral(nt)) * M ...
                + structure.FDMatrix.Ay2 .* exp(-structure.rpbcFactorY .* pulse.gIntegral(nt)) * M ...
                + structure.FDMatrix.Az1 .* exp( structure.rpbcFactorZ .* pulse.gIntegral(nt)) * M ...
                + structure.FDMatrix.Az2 .* exp(-structure.rpbcFactorZ .* pulse.gIntegral(nt)) * M ; 
        else                         % no periodic boundary condition
            M = structure.FDMatrix.A * M ; 
        end

        % T2 relaxation
        if structure.useT2 == 'y'     % T2 needs to be considered
            M = structure.T2relaxation .* M ; 
        end    
    end

       
end


