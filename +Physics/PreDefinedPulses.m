function pulse = PreDefinedPulses(AcqMethod) 
%% This function loads pre-defined pulse sequence that have been published previously. 
% There are currently four options
% 1. ADCSpectrum: multiple diffusion times ranging from 1.25 ms to 1000ms with the same b value
%       Reference: Xu et al. 
% 2. The IMPULSED protocol for animal imaging with Gmax<=360mT/m on each gradietn axis. Usually all gradient axes are turned on 
%       Reference: Jiang et al. Magn Reson Med. 2017;78(1):156-64. Epub 2016/08/09. doi: 10.1002/mrm.26356. PubMed PMID: 27495144; 
% 3. The IMPULSED protocol for human imaging with Gmax<=80mT/m for all gradients on. 
%       Reference: Xu et al. Magn Reson Med. 2020;83(6):2002-14. Epub 20191125. doi: 10.1002/mrm.28056. PubMed PMID: 31765494; PMCID: PMC7047520.
% 4. The VERDICT protocol for human imaging particularly for prostate imaging. 
%       Reference: Johnston et al. Radiology. 2019;291(2):391-7. Epub 20190402. doi: 10.1148/radiol.2019181749. PubMed PMID: 30938627; PMCID: PMC6493214.
%
% ---------------------------------------------------------------------------------------------------------------------
% NOTE
%       1. Although it is flexible to specify any sequence parameters, it is recommend to set the first acquisition as b=eps as the T2-weighted b=0 signal.
% 
% ---------------------------------------------------------------------------------------------------------------------
% Author: Junzhong (JZ) Xu, Dec 25, 2022
% 
% ---------------------------------------------------------------------------------------------------------------------


%%create pulse using the mati package
    switch lower(AcqMethod)
        case 'adc_spectrum'      %------------- ADC spectrum with multiple diffusion times but with the same b value -------------------
            Nacq = 30 ;          % total number of acquisition points
            pulse = mati.DiffusionPulseSequence(Nacq,...
                'TE',               110,...                 % echo time [ms]
                'delta',            [10,15,20, 25, 30, 35, 40, 20, 25, 30, 35, 40, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,12,12,12,12],...                  % gradient duration [ms]
                'Delta',            [40,40,40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 14, 19, 24, 34, 54, 74, 104, 154, 204, 254, 304,404,504,604,704,804,904,1004],...                  % separation of two gradients [ms]
                'b',                  [ones(1,30)*0.2], ...      % b value [ms/μm^2]
                'n',                  [2,2,2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0],...             % number of oscillating cycles
                'shape',         ["tcos","tcos", "tcos","tcos", "tcos","tcos", "tcos", "tcos", "tcos", "tcos", "tcos", "tcos", "tpgse","tpgse", "tpgse","tpgse", "tpgse","tpgse", "tpgse","tpgse", "tpgse", "tpgse", "tpgse","tpgse", "tpgse", "tpgse","tpgse", "tpgse", "tpgse","tpgse"],...                 % gradient waveform shape
                'gdir',             [1 0 0],...              % gradient directions. It should be a Nx3 matrix
                'trise',            0.1) ;                  % gradient rise time [ms]
        case 'impulsed_animal'       % ------------ the IMPULSED protocol for animal imaging with Gmax<=360 mT/m on each gradient axis ------------------
            Nacq = 16 ;          % total number of acquisition points
            pulse = mati.DiffusionPulseSequence(Nacq,...
                'TE',               [67] ,...                 % echo time [ms]
                'delta',            [ones(1,4)*20, ones(1,4)*20, ones(1,4)*20, ones(1,4)*3],...                  % gradient duration [ms]
                'Delta',            [ones(1,4)*25, ones(1,4)*25, ones(1,4)*25, ones(1,4)*46],...                  % separation of two gradients [ms]
                'b',                  [[0.15, 0.3,0.45,0.6], [0.33, 0.66, 1.0, 1.32], [0.375, 0.750, 1.125, 1.50], [0.375, 0.750, 1.125, 1.50]], ...      % b value [ms/μm^2]
                'n',                  [ones(1,4)*3, ones(1,4)*2, ones(1,4)*1, zeros(1,4)],...             % number of oscillating cycles
                'shape',         [repmat("cos",[1,4]), repmat("cos",[1,4]), repmat("cos",[1,4]), repmat("tpgse",[1,4])],...                 % gradient waveform shape
                'gdir',             [1 0 0],...              % gradient directions. It should be a Nx3 matrix
                'trise',            0.9) ;                  % gradient rise time [ms]
        case 'impulsed_human'       % ------------ the IMPULSED protocol for human imaging with Gmax<=80 mT/m of all gradients -------------------
            Nacq = 13 ;          % total number of acquisition points
            pulse = mati.DiffusionPulseSequence(Nacq,...
                'TE',               110,...                 % echo time [ms]
                'delta',            [ones(1,3)*40.9, ones(1,4)*40.9, ones(1,6)*12],...                  % gradient duration [ms]
                'Delta',            [ones(1,3)*51.4, ones(1,4)*51.4, ones(1,6)*74],...                  % separation of two gradients [ms]
                'b',                  [[0.1, 0.2, 0.3], 0.25:0.25:1, [0.25:0.25:1, 1.4, 1.8]], ...      % b value [ms/μm^2]
                'n',                  [ones(1,3)*2, ones(1,4), zeros(1,6)],...             % number of oscillating cycles
                'shape',         [repmat("tcos",[1,3]), repmat("tcos",[1,4]), repmat("tpgse",[1,6])],...                 % gradient waveform shape
                'gdir',             [1 0 0],...              % gradient directions. It should be a Nx3 matrix
                'trise',            0.9) ;                  % gradient rise time [ms]
        case 'verdict'
            Nacq = 5 ;          % total number of acquisition points
            pulse = mati.DiffusionPulseSequence(Nacq,...
                'TE',               [50, 65, 90, 71, 80],...                 % echo time [ms]
                'delta',            [3.9, 11.4, 23.9, 14.4, 18.9],...                  % separation of two gradients [ms]
                'Delta',            [23.8, 31.3, 43.8, 34.3 38.8],...                  % gradient duration [ms]
                'b',                  [0.09, 0.5, 1.5, 2.0, 3.0], ...      % b value [ms/μm^2]
                'n',                  [0],...             % number of oscillating cycles
                'shape',         ["tpgse"],...                 % gradient waveform shape
                'gdir',             [1 0 0],...              % gradient directions. It should be a Nx3 matrix
                'trise',            0.9) ;                  % gradient rise time [ms]
        case 'ssift'
            Nacq = 2 ;          % total number of acquisition points
            pulse = mati.DiffusionPulseSequence(Nacq,...
                'TE',               [110],...                 % echo time [ms]
                'delta',            [40.9, 12],...                  % separation of two gradients [ms]
                'Delta',            [51.4, 74],...                  % gradient duration [ms]
                'b',                  [1, 1], ...      % b value [ms/μm^2]
                'n',                  [1, 0],...             % number of oscillating cycles
                'shape',         ["tcos", "tpgse"],...                 % gradient waveform shape
                'gdir',             [1 0 0],...              % gradient directions. It should be a Nx3 matrix
                'trise',            0.9) ;                  % gradient rise time [ms]
        otherwise
            error('The input nMethod should be 1 (ADCSpectrum), 2 (IMPULSEDAnimal), 3(IMPULSEDHuman), or 4(VERDICT)')
    end
end
