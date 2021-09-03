%% Example of calculating ADC spectra of one restricted compartment
% This example shows how to calculate frequency-dependent OGSE ADC spectra of one restricted compartment based on analytical equations. 
% The restricted compartments can be parallel planes, cylinders, spheres, or spherical shells (hollow spheres). 
% NOTE: diffusion is perfectly restricted without surface relaxation or permeability. 
% There are two independent microstructural parameters: 
%
% # compartment size $d$: distance between planes, diameters of cylinders or spheres. NOTE: spherical shells need two diameters din and dout. 
% # intracompartment diffusion coefficient: $D_{in}$
%
% 
% *Reference* 
% 
% The main reference that should be cited when using the code in this script is
% 
% # Xu J, et al. Quantitative characterization of tissue microstructure with temporal diffusion spectroscopy. J Magn Reson. 2009;200(2):189-97. PubMed PMID: 19616979.
% # Xu J, et al. Magnetic resonance imaging of mean cell size in human breast tumors. Magn Reson Med. 2020;83(6):2002-14. PubMed PMID: 31765494.
% 
% *Comments or questions?* 
% 
% Please send your comments or questions to Junzhong (JZ) Xu (junzhong.xu@vanderbilt.edu)


%% Preliminary
clear ; clear obj ; 

%% Generate a PulseSequence object
% Set pulse parameters 
delta = 20 ;                            % each gradient duration [ms]
fs = [25:25:1000]*1e-3 ;        % the range of gradient frequencies [kHz]
ns = floor(delta .* fs) ;           % the range of number of oscillating cycles 
Nacq = length(ns) ;               % total number of data acquisition points

% Create the object
pulse = mati.DiffusionPulseSequence(sum(Nacq) , ...
                    'delta',            delta , ...         % gradient duration [ms]
                    'Delta' ,          delta+5 , ...     % gradient separation [ms]
                    'shape',         "cos" , ...        % gradient waveform shape, a single string or a string array
                    'b' ,               1 , ...               % b value [ms/um^2]
                    'n' ,               ns ...               % number of gradient oscillating cycles
            ) ; 

       
%%  Example#1: ADC spectrum of restricted diffusion between parallel planes with different distances

% Generate a structure object
structure.modelName = '1compt' ; structure.geometry = 'plane' ; 
plane_model = mati.IMPULSED(structure, pulse) ; 

% Set microstructual parameters that are of interest
d = [2:2:10] ; Din = [1.56] ; 

% Calculate restricted dMRI signals based on analytical equations
signal_plane = plane_model.FcnSignal({d, Din}, plane_model) ; 

% Calculate ADC based on signals and b values
ADC_plane = -log(signal_plane) ./ pulse.b ; 

% Plot ADC spectra of restriction between parallel planes with different distances
figure(1) ; clf ; 
plot(pulse.f*1000, ADC_plane, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('d =  2\mum','d =  4\mum','d =  6\mum','d =  8\mum','d = 10\mum','Location','Southeast') ; title('Plane, Varying d')



%%  Example#2: ADC spectrum of restricted diffusion perpendicular to cylinders with different intra-cylinder diffusion coefficients

% Generate a structure object
structure.modelName = '1compt' ; structure.geometry = 'cylinder' ; 
cylinder_model = mati.IMPULSED(structure, pulse) ; 

% Set microstructural parameters of interest
d = [5] ; Din = [1:0.2:2] ; [d,Din] = meshgrid(d,Din) ; sim_parms = [d(:)'; Din(:)'] ;     

% Calculate restricted dMRI signals based on analytical equations
signal_cylinder = cylinder_model.FcnSignal(sim_parms, cylinder_model) ; 

% Calculate ADC values
ADC_cylinder = -log(signal_cylinder) ./ pulse.b ; 

% Plot ADC spectra of restricted diffusion perpendicular to cylinders with different intra-cylinder diffusion coefficients
figure(2) ; clf ; 
plot(pulse.f*1000, ADC_cylinder, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('1\mum^2/ms','1.2\mum^2/ms','1.4\mum^2/ms','1.6\mum^2/ms','1.8\mum^2/ms','Location','Southeast'); title('Cylinder, Varying Din')


%%  Example#3: ADC spectra of restricted diffusion inside spheres with different diameters

% Generate a structure object
structure.modelName = '1compt' ; structure.geometry = 'sphere' ; 
sphere_model = mati.IMPULSED(structure, pulse) ; 

% Set microstructural parameters of interest
d = [5:5:25] ; Din = [1.56] ; [d,Din] = meshgrid(d,Din) ; sim_parms = [d(:)'; Din(:)'] ;     

% Calculate restricted dMRI signals based on analytical equations
signal_sphere = sphere_model.FcnSignal(sim_parms, sphere_model) ; 

% Calculate ADC values
ADC_sphere = -log(signal_sphere) ./ pulse.b ; 

% Plot ADC spectra of restricted diffusion inside spheres with different diameters 
figure(3) ; clf ; 
plot(pulse.f*1e3, ADC_sphere, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('d =  5\mum','d = 10\mum','d = 15\mum','d = 20\mum','d = 25\mum','Location','Southeast') ; title('sphere, Varying d')


%% Example#4: ADC spectrum of restricted diffusion inside spherical shells (hollow spheres) with different din values

% Generate a structure object
structure.modelName = '1comptHollow' ; structure.geometry = 'hollowSphere' ; 
shell_model = mati.IMPULSED(structure, pulse) ; 

% Set microstructural parameters of interest 
din = [1:3:15] ; dout = [16] ; Din = [1.56] ; [din, dout,Din] = meshgrid(din, dout,Din) ; sim_parms = [din(:)'; dout(:)'; Din(:)'] ; 
signal_shell = shell_model.FcnSignal(sim_parms, shell_model) ; 

% Calculate ADC values
ADC_shell = -log(signal_shell) ./ pulse.b ; 

% Plot ADC spectra of restricted diffusion inside spherical shells with different din values
figure(4) ; clf ; 
plot(pulse.f*1000, ADC_shell, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('din = 1\mum','din = 4\mum','din = 7\mum','din = 10\mum','din = 13\mum','Location','Southeast') ; title('Hollow Sphere, dout=15\mum, Varying din')



