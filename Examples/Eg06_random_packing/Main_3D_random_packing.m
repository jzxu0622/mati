%% Create 3D structure of randomly packed spheres
% the main script code is based on https://github.com/VasiliBaranov/packing-generation. 
% Please copy this function to packing-generation-master\Run\matlab to run
% citation: 
% 1. Baranau and Tallarek (2014) Random-close packing limits for monodisperse and polydisperse hard spheres, doi:10.1039/C3SM52959B. 
% 2. Or, alternatively, Baranau et al. (2013) Pore-size entropy of random hard-sphere packings; doi:10.1039/C3SM27374A.
% 
% -----------------------------------------------------------------------------------------------------------------------------------------

clear
coreFolder = '.' ; 

%% create cell diameter distribution
% !!! adjust the following parameters to decide the tissue size, # of spheres, and the cell size distribution
Nx = 128 ;      Ny = Nx ; Nz = Nx ;     % size of total 3D structure in each direction
Nd = 2700 ;         % # of spheres
dmean = 10 ;        % mean sphere diameter
dsigma = 2 ;         % std of sphere diameter

d = round(abs(randn([Nd 1])*dsigma+dmean)) ; 
if d<1, d=1 ; end % dmin =1 for discretization

vin0 = sum(4/3*pi*(d/2).^3) / (Nx*Ny*Nz) ; 

% visualization
figure(1) ; clf ; 
histogram(d)
title(sprintf('vin  = %.02f', vin0)) ; 

%  output to diameters.txt
filename = fopen(fullfile(coreFolder,'diameters.txt'),'w') ; 
fprintf(filename, '%d\n', d) ; 
fclose(filename) ; 

%% update generation.conf
filename = fopen('generation.conf','w') ; 
fprintf(filename, 'Particles count: %d\n', Nd) ; 
fprintf(filename, 'Packing size: %.1f %.1f %.1f\n', Nx, Ny, Nz) ; 
fprintf(filename, 'Generation start: 1\n') ;  
fprintf(filename, 'Seed: %d\n', randi([1 1000],1,1)) ; 
fprintf(filename, 'Steps to write: 1000\n') ; 
fprintf(filename, 'Boundaries mode: 1\n') ; 
fprintf(filename, 'Contraction rate: 1.328910e-005\n') ; 
fprintf(filename, '1. boundaries mode: 1 - bulk; 2 - ellipse (inscribed in XYZ box, Z is length of an ellipse); 3 - rectangle\n') ; 
fprintf(filename, '2. generationMode = 1 (Poisson, R) or 2 (Poisson in cells, S)\n') ; 
fclose(filename) ; 

%% Run VasiliBaranov's code
disp('-------------- start packing -------------------')
system([fullfile(coreFolder, 'PackingGeneration.exe'), ' -fba | tee log_fba.txt']) ; 
system([fullfile(coreFolder, 'PackingGeneration.exe'), ' -ls | tee log_ls.txt']) ; 
system([fullfile(coreFolder, 'PackingGeneration.exe'), ' -lsgd | tee log_lsgd.txt']) ; 
disp('-------------- end of packing -------------------')


%% Read packing results
packingSerializer = PackingSerializer();

% --------------- Reading the config -------------------
configFilePath = fullfile(coreFolder, 'generation.conf');
config = packingSerializer.ReadConfig(configFilePath);

% -------------------- Reading the packing -----------------------
packingFilePath = fullfile(coreFolder, 'packing_init.xyzd');
[packing, floatDataType] = packingSerializer.ReadPacking(packingFilePath, config.ParticlesCount);

% One can also just call packing = packingSerializer.ReadPacking(packingFilePath);
% config.ParticlesCount helps to determine the real precision with which particles were stored (single or double), 
% otherwise it is assumed to be double. The C++ code on github also uses double precision.

%%%%%%%%%%%%%%
% IMPORTANT!!!
% Final diameters will usually be smaller than the ones you specify in the original packing 
% (they are linearly scaled at the beginning of generation by a small value (scaling factor) 
% and this value increases during generation, but usually not until 1). 
% The final packing (packing.xyzd) will store correct particle centers and the ORIGINAL diameters 
% (not the ones scaled with the final scaling factor)for historical reasons. 
% To get the final scaling factor and correct final diameters, 
% you have to use the Final Porosity field from the packing.nfo. 
% The equation for the final scaling factor is 
% FinalDensity = 1 - FinalPorosity = N * pi/6 * (D_original * FinalScalingFactor)^3 / (Lx * Ly * Lz) = 
% TheoreticalDensity*FinalScalingFactor^3 = (1-TheoreticalPorosity)*FinalScalingFactor^3. 
% Thus, FinalScalingFactor = ((1 - FinalPorosity) / (1-TheoreticalPorosity))^(1/3).
% The equation is the same for polydisperse packings as well.
% Theoretical porosity is also written in packing.nfo.
% When you scale diameters with the final scaling factor, you don't need to scale the box size or particle centers.

% Reading packing.nfo
infoFilePath = fullfile(coreFolder, 'packing.nfo');
packingInfo = packingSerializer.ReadPackingInfo(infoFilePath);

% Scaling the diameters. You can of course scale the coordinates and the box size instead (with 1/finalScalingFactor).
finalScalingFactor = ((1 - packingInfo.FinalPorosity) / (1 - packingInfo.TheoreticalPorosity))^(1/3);
packing.ParticleDiameters = packing.ParticleDiameters * finalScalingFactor;

% Now you can really use packing.ParticleCoordinates and packing.ParticleDiameters

packing.BoxDimensions = config.BoxDimensions; % just for the sake of completeness

% If you want to save the updated packing, use smth like
% packingSerializer.WritePacking(packingFilePath, packing, floatDataType);
% % packingSerializer.WriteConfig(configFilePath, config); % needed only if you scale particle centers and the box size instead of diameters
% You may also want to update packing.nfo and set TheoreticalPorosity to FinalPorosity to avoid scaling the packing once again 
% the next time you run this script. Or use another name for the scaled packing

%% vin
center = round([packing.ParticleCoordinates, packing.ParticleDiameters/2, (1:length(packing.ParticleDiameters))']) ; 
Nx = packing.BoxDimensions(1) ; Ny = packing.BoxDimensions(2) ; Nz = packing.BoxDimensions(3) ; 
r = center(:,4) ; 
% vin  = sum(4/3*pi*(r).^3) / prod(packing.BoxDimensions)

% enlarge spheres a bit so enhance intra volume fraction. This parameter needs to be manually adjusted to form a model tissue that is wanted. 
r_enlarge = 1.1 ; 

% %% save matlab data
% save(fullfile(folder,'Result_packing.mat')) ; 

%% output tissue file
tissue = zeros(Nx,Ny,Nz); 
for num=1:size(center,1)
    for i=1:Nx
        for j = 1:Ny
            for k = 1:Nz
                if (i-center(num,1))^2 + (j-center(num,2))^2 + (k-center(num,3))^2 <= (center(num,4)*r_enlarge)^2 
                    tissue(i,j,k) = center(num,5) ; 
                end
            end
        end
    end
end


%% calculate the volumn fractions 
vi = sum(tissue(:) >= 1) / (Nx*Ny*Nz) ;
ve = sum(tissue(:) == 0) / (Nx*Ny*Nz) ;

fprintf(' # of spheres: %d\n', size(center,1)) ;  
fprintf('r_enlarge = %f\n', r_enlarge) ; 
fprintf('vi = %.02f, ve = %.02f\n', vi, ve)   

disp('----------------- Press any key if the vi is ok. Otherwise, Ctrl+C ------------------')
pause
disp('OK! writing structure file ...')

%% save 3D tissues to 2D
[Nx, Ny, Nz] = size(tissue) ; 
tissue2d = zeros([Nx*Nz Ny]) ; 
for nz=1:Nz
    tissue2d(((nz-1)*Nx+1):(nz*Nx),:) = tissue(:,:,nz) ; 
end
filename = sprintf('3Dvin=%02d_%d_%d_%d.dat', round(vi*100),Nx, Ny, Nz) ; 
writematrix(tissue2d, filename,'Delimiter',' ')

% output tissue file
save(sprintf('3Dvin=%02d_%d_%d_%d.mat', round(vi*100), Nx, Ny, Nz)) 


%% visualization of 3D packed cells
disp(' data files written. Now for visualization') 
figure(2) ; clf ; hold on ; 
for n=1:size(center,1)
    [x,y,z]=ellipsoid(center(n,1),center(n,2),center(n,3),center(n,4),center(n,4),center(n,4)*r_enlarge); 
    hsph = surf(x,y,z) ; 
    hsph.EdgeColor = 'none' ;     hsph.FaceColor = 'g' ; 
end
daspect([1 1 1])
axis equal
axis tight
camlight 
lighting phong
axis off
view(126,17)
material shiny

Xcorner = [0,   0,   Nx,   Nx, 0, 0, 0,   Nx,   Nx, 0, ]  ; 
Ycorner = [0,   Ny, Ny,  0, 0, 0, Ny, Ny,  0, 0, ] ; 
Zcorner = [0,   0,   0,     0, 0, Nz, Nz, Nz, Nz, Nz]  ; 
plot3(Xcorner, Ycorner, Zcorner,'k', 'LineWidth',2) ; 
plot3([Nx Nx],[0 0],[0 Nz], 'k','Linewidth',2) 
plot3([0 0],[Ny Ny],[0 Nz], 'k','Linewidth',2) 
plot3([Nx Nx],[Ny Ny],[0 Nz], 'k','Linewidth',2) 

%% End
disp('-------------- End of Main_3D_random_packing.m -----------------')
