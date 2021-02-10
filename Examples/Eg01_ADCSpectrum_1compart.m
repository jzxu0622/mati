%% Example of calculating ADC spectra of one restricted compartments
% 
% -------------------------------------------------------------------------------------------------------------------------

%% Preliminary
clear ; clear obj ; 

%% Generate pulse object
delta = 20 ; % [ms]
fs = [25:25:1000]*1e-3 ; % [kHz]
ns = floor(delta .* fs) ; 
N = length(ns) ;  
pulse = mig.DiffusionPulseSequence(sum(N) , ...
                    'delta',           delta , ...         % gradient duration [ms]
                    'Delta' ,         delta+5 , ...     % gradient separation [ms]
                    'shape',         "cos" , ...        % gradient shape, a single string or a string array
                    'b' ,                1 , ...               % b value [ms/um^2]
                    'n' ,                 ns ...               % number of gradient oscillating cycles
            ) ; 

       
%%  ADC spectrum of planes
% generate structure object
structure.modelName = '1compt' ; structure.geometry = 'plane' ; 
plane_model = mig.IMPULSED(structure, pulse) ; 

% plane distance array
d = [2:2:10] ; Din = [1.56] ; 
% calculate restricted DWI signal
signal_plane = plane_model.FcnSignal({d, Din}, plane_model) ; 
% signal_plane = plane_model.FcnSignal([2 4 6 8 10; 1 2 1 2 1], plane_model) ; 
ADC_plane = -log(signal_plane) ./ pulse.b ; 

figure(1) ; clf ; 
plot(pulse.f*1000, ADC_plane, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('d =  2\mum','d =  4\mum','d =  6\mum','d =  8\mum','d = 10\mum','Location','Southeast') ; title('Plane, Varying d')

%%  ADC spectrum of cylinders
% generate structure object
structure.modelName = '1compt' ; structure.geometry = 'cylinder' ; 
cylinder_model = mig.IMPULSED(structure, pulse) ; 

% plane distance array
d = [5] ; Din = [1:0.2:2] ; [d,Din] = meshgrid(d,Din) ; sim_parms = [d(:)'; Din(:)'] ;     
% calculate restricted DWI signal
signal_cylinder = cylinder_model.FcnSignal(sim_parms, cylinder_model) ; 
ADC_cylinder = -log(signal_cylinder) ./ pulse.b ; 

figure(2) ; clf ; 
plot(pulse.f*1000, ADC_cylinder, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('1\mum^2/ms','1.2\mum^2/ms','1.4\mum^2/ms','1.6\mum^2/ms','1.8\mum^2/ms','Location','Southeast'); title('Cylinder, Varying Din')

%%  ADC spectrum of sphere
% generate structure object
structure.modelName = '1compt' ; structure.geometry = 'sphere' ; 
sphere_model = mig.IMPULSED(structure, pulse) ; 

% plane distance array
d = [5:5:25] ; Din = [1.56] ; [d,Din] = meshgrid(d,Din) ; sim_parms = [d(:)'; Din(:)'] ;     
% calculate restricted DWI signal
signal_sphere = sphere_model.FcnSignal(sim_parms, sphere_model) ; 
ADC_sphere = -log(signal_sphere) ./ pulse.b ; 

figure(3) ; clf ; 
plot(pulse.f*1e3, ADC_sphere, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('d =  5\mum','d = 10\mum','d = 15\mum','d = 20\mum','d = 25\mum','Location','Southeast') ; title('sphere, Varying d')


%% ADC spectrum of hollow sphere
structure.modelName = '1comptHollow' ; structure.geometry = 'hollowSphere' ; 
shell_model = mig.IMPULSED(structure, pulse) ; 

din = [1:3:15] ; dout = [16] ; Din = [1.56] ; [din, dout,Din] = meshgrid(din, dout,Din) ; sim_parms = [din(:)'; dout(:)'; Din(:)'] ;     
signal_shell = shell_model.FcnSignal(sim_parms, shell_model) ; 
ADC_shell = -log(signal_shell) ./ pulse.b ; 

figure(4) ; clf ; 
plot(pulse.f*1000, ADC_shell, 'o-') ; xlabel('f [Hz]') ; ylabel('ADC [\mum^2/ms]') ; 
legend('din = 1\mum','din = 4\mum','din = 7\mum','din = 10\mum','din = 13\mum','Location','Southeast') ; title('Hollow Sphere, dout=15\mum, Varying din')

%% End

