%% Evaluate the short pulse approximaiton by comparing short vs finite gradients
% 
% -------------------------------------------------------------------------------------------------------------------------


%% Preliminary
clear ; clear obj ; 

%% Generate pulse object
Deltas = [4:20] ; 
deltas = Deltas/2 ; % [ms]
N = length(Deltas) ;  
pulse = mati.DiffusionPulseSequence(sum(N) , ...
                    'delta',           deltas , ...         % gradient duration [ms]
                    'Delta' ,         Deltas , ...     % gradient separation [ms]
                    'shape',         "pgse" , ...        % gradient shape, a single string or a string array
                    'b' ,                1 ...               % b value [ms/um^2]
            ) ; 

       
%%  ADC of long pulses
% generate structure object
structure.modelName = '1compt' ; structure.geometry = 'plane' ; 
sphere_model = mati.IMPULSED(structure, pulse) ; 

% plane distance array
d = [10:2:16] ; Din = [1.56] ; 
% calculate restricted DWI signal
signal_signal = sphere_model.FcnSignal({d, Din}, sphere_model) ; 
% signal_plane = plane_model.FcnSignal([2 4 6 8 10; 1 2 1 2 1], plane_model) ; 
ADC_plane = -log(signal_signal) ./ pulse.b ; 


%% Generate pulse object
deltas = 2 ; % [ms]
N = length(Deltas) ;  
pulse2 = mati.DiffusionPulseSequence(sum(N) , ...
                    'delta',           deltas , ...         % gradient duration [ms]
                    'Delta' ,         Deltas , ...     % gradient separation [ms]
                    'shape',         "pgse" , ...        % gradient shape, a single string or a string array
                    'b' ,                1 ...               % b value [ms/um^2]
            ) ; 

       
%%  ADC of long pulses
% generate structure object
structure.modelName = '1compt' ; structure.geometry = 'sphere' ; 
sphere_model2 = mati.IMPULSED(structure, pulse2) ; 

% calculate restricted DWI signal
signal_sphere2 = sphere_model2.FcnSignal({d, Din}, sphere_model2) ; 
% signal_plane = plane_model.FcnSignal([2 4 6 8 10; 1 2 1 2 1], plane_model) ; 
ADC_plane2 = -log(signal_sphere2) ./ pulse.b ; 


%%
figure(1) ; clf ; 
plot(pulse.Delta, signal_signal,'o') ; hold on; 
plot(pulse.Delta,signal_sphere2, '-') ; xlabel('Delta [ms]') ; ylabel('signal ') ; 
legend('d =  10\mum','d = 12\mum','d =  14\mum','d =  16\mum','Location','Southeast') ; title('short vs finite gradient duration')

%% End

