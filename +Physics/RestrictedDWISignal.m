function E = RestrictedDWISignal(parms, model)
%% Calculates analytical diffusion MRI signal attenuation of restricted diffusion 
% The geometry must be a member of {plane, cylinder, sphere, or spherical shell}.
% The separation of pulse sequence and microtructure parameters is based on "Janez Stepisnik, Physica tp, 183, 343-350 (1993)"
% The signal equations of sin or cos and microstructure expressions of spherical shell are based on " Xu J et al. Quantitative characterization of tissue microstructure with temporal diffusion spectroscopy. J Magn Reson 2009;200(2):189-197."
% The signal equations of tpgse or tcos shapes are based on "Xu J et al. Magnetic resonance imaging of mean cell size in human breast tumors. Magn Reson Med. 2020 Jun;83(6):2002-2014."
% 
% -------------------------------------------------------------------------------------------------------------------------
% INPUTS: 
%       model:
%           model.pulse: pulse sequence parameters. The elements such as f, T, delta, tau, w should be scalars or matrices with the same dimensions
%               model.pulse.shape: pgse, tpgse, sin, cos, tcos, tsin
%           model.structure: microstructural parameters: geometry can be plane, cylinder, sphere, or spherical shell
%       parms: 
%           for {'plane','cylinder','sphere'}, parms is a 2 x N marix [d; D]
%           for {'hollowSphere', 'sphericalShell'}, parms is a 3 x N matrix [din; dout; D]
% OUTPUTS:
%       E:   signal attenuation with the same dimension of pulse elements
% -------------------------------------------------------------------------------------------------------------------------
% NOTES:
%        1. For plane: d is the distance between planes
%        2. For cylinder or sphere: d is the diameter
%        3. All inputs are d (diameter)
%        4. The model.pulse can take any combinations of pulse sequences (pgse, tpgse, sin, cos, tcos, ...)
%        5. model.pulse.shape == cos is NOT apodised but an ideal cosine shape !!!!!!! It provides little bias for fitting size 
%        6. xstep in the function get_hollow_sphere_roots(a,b) should be small enough to count for a small din (inner diameter)
% -------------------------------------------------------------------------------------------------------------------------
% Author: Junzhong Xu, Dec 08, 2020
% 
% -------------------------------------------------------------------------------------------------------------------------

    %% check dimensions: d and D should be either scalars or row vectors of the same size
    structure = model.structure ; pulse = model.pulse ; 
%     pulse
%     structure
    switch structure.geometry
        case {'plane','cylinder','sphere'}
            d = parms(1,:) ; D = parms(2,:) ; 
            % check inputs
                if ~isvector(d) || ~isvector(D), error('%s: The input \{d, D\} should be vectors of the same size', mfilename) ; end
                if ~isscalar(d) && iscolumn(d), d = d' ; end
                if ~isscalar(D) && iscolumn(D), D = D' ; end
                if isscalar(d) && ~isscalar(D), d = repmat(d, [1 length(D)]) ; end
                if isscalar(D) && ~isscalar(d), D = repmat(D, [1 length(d)]) ; end
                if any(length(d)~=length(D)), error('%s: both d and D should be vectors of the same size', mfilename) ; end
            % check inputs
                if any(d<=0), error('%s: The diameter d should be a positive number', mfilename); end
                if any(D<=0), error('%s: D should be a positive number ans smaller than the diffusion coefficient of free water', mfilename); end
%                 if any(D>3.07*1.1), warning('%s: D is larger than 110% of 3.07um^2/ms, the diffusion coefficient of free water at 37C', mfilename); end
            % convert to R (radius)
                R = d/2 ; 
            % load Bk and ak
                [Bk, ak] = load_Bk_ak(R, structure) ; 
        case {'hollowSphere', 'sphericalShell'}
            din = parms(1,:) ; dout = parms(2,:) ; D = parms(3,:) ; 
            % check inputs
                if any(din>=dout), error('din (inner diameter) should be smaller than dout (outer diameter)') ; end
            % convert d (diameter) to R (radius)
                Rin = din/2 ; Rout = dout/2 ; 
            % load Bk and ak
                [Bk, ak] = load_Bk_ak_hollow(Rin, Rout) ; 
        otherwise
            error('%s: The input model.structure.geometry should be \{''plane'',''cylinder'',''sphere'',''hollow_sphere'',''spherical_shell''\}',mfilename) ; 
    end

    
    %% preliminary
    % convert diameter to radius. All followings use radius
    for n=1:length(pulse.shape)
        if ~ismember(pulse.shape(n), ["pgse"; "tpgse"; "sin"; "cos"; "tcos"; "tsin"] ) 
            error('%s: The #%d shape (%s) is not supported. A shape should be one of [pgse, tpgse, sin, cos, tcos, tsin]', mfilename, n, pulse.shape(n)); 
        end
    end
    
    % get indices of each gradient shape
    Ind.pgse = (pulse.shape == 'pgse') ; 
    Ind.tpgse = (pulse.shape == 'tpgse') ; 
    Ind.sin = (pulse.shape == 'sin') ; 
    Ind.cos = (pulse.shape == 'cos') ; 
    Ind.tcos1 = (pulse.shape == 'tcos' & pulse.n==1) ; 
    Ind.tcos2 = (pulse.shape == 'tcos' & pulse.n==2) ; 
    Ind.tcos3 = (pulse.shape == 'tcos' & pulse.n==3) ; 
    Ind.tsin1 = (pulse.shape == 'tsin' & pulse.n==1) ; 
    Ind.tsin2 = (pulse.shape == 'tsin' & pulse.n==2) ; 

    %% calculate structural parameters
    [~, NR, Nk] = size(Bk) ; 
    Bk = repmat(Bk, [pulse.Nacq, 1, 1]) ; ak = repmat(ak, [pulse.Nacq, 1,1]) ; 
    D = repmat(D, [pulse.Nacq, 1, Nk]) ; 
    % reshape pulse
    delta = pulse.delta ; Delta = pulse.Delta  ;             
    delta = repmat(delta,[1 NR, Nk]) ; Delta = repmat(Delta,[1 NR, Nk]) ; 
    
    %% calculate attenuation
    % First, pgse for all as a place holder
    Etmp = 1./ak.^2./D.^2.*(ak.*D.*delta - 1 + exp(-ak.*D.*delta) + exp(-ak.*D.*Delta)- 1/2*exp(-ak.*D.*(Delta+delta)) - 1/2*exp(-ak.*D.*(Delta-delta))) ; 
    % tpgse
    if any(Ind.tpgse)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tpgse = 1/4*(12*exp(-(ak.*D).*tp) - 24*exp(-(ak.*D).*trise) - 24*exp(-Delta.*(ak.*D)) - 24*(ak.*D).*trise + 12*exp(-(ak.*D).*(Delta-trise)) + 12*exp(-(ak.*D).*(Delta+trise)) - 6*exp(-(ak.*D).*(Delta-tp)) - 6*exp(-(ak.*D).*(Delta+tp)) - 24*exp(-(ak.*D).*(trise+tp)) + 12*exp(-(ak.*D).*(2*trise+tp)) + 8*(ak.*D).^3.*trise.^3 + 12*(ak.*D).^3.*trise.^2.*tp + 12*exp(-(ak.*D).*(Delta-trise-tp)) + 12*exp(-(ak.*D).*(Delta+trise+tp)) - 6*exp(-(ak.*D).*(Delta-2*trise-tp)) - 6*exp(-(ak.*D).*(Delta+2*trise+tp)) + 24)./(3*(ak.*D).^4.*trise.^2) ;
        Etmp(Ind.tpgse,:,:) = Etmp_tpgse(Ind.tpgse,:,:) ; 
    end
    % sin
    if any(Ind.sin)
        w=pulse.w ; w = repmat(w,[1 NR, Nk]) ; 
        Etmp_sin = w.^2 ./(ak.^2.*D.^2 + w.^2).^2 .* (ak.*D.*(ak.^2.*D.^2 + w.^2).*delta./2./w.^2 +1 - exp(-ak.*D.*delta) - exp(-ak.*D.*Delta)  +0.5.*exp(-ak.*D.*(delta+Delta)) +0.5.*exp(-ak.*D.*(Delta-delta))  ) ;  
        Etmp(Ind.sin,:,:) = Etmp_sin(Ind.sin,:,:) ; 
    end
    % cos
    if any(Ind.cos)
        w=pulse.w ; w = repmat(w,[1 NR, Nk]) ; 
        Etmp_cos = ak.^2 .*D.^2 ./(ak.^2.*D.^2 + w.^2).^2 .* ((ak.^2.*D.^2 + w.^2).*delta./2./ak./D -1 + exp(-ak.*D.*delta) + exp(-ak.*D.*Delta)   -0.5.*exp(-ak.*D.*(delta+Delta)) -0.5.*exp(-ak.*D.*(Delta-delta))  ) ; 
        Etmp(Ind.cos,:,:) = Etmp_cos(Ind.cos,:,:) ; 
    end
    % tcos and n=1
    if any(Ind.tcos1)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tcos1 = 1/4*(8*exp(-(ak.*D).*tp) - 8*exp(-(ak.*D).*trise) - 8.*exp(-2.*(ak.*D).*trise) - 16.*exp(-(ak.*D).*Delta) - 24.*(ak.*D).*trise + 4.*exp(-(ak.*D).*(Delta-trise)) + 4.*exp(-(ak.*D).*(Delta+trise)) + 4.*exp(-(ak.*D).*(Delta+2.*trise)) + 4.*exp(-(ak.*D).*(Delta-2.*trise)) - 4.*exp(-(ak.*D).*(Delta-tp)) - 4.*exp(-(ak.*D).*(Delta+tp)) - 8.*exp(-(ak.*D).*(trise+tp)) + 4.*exp(-(ak.*D).*(trise+2.*tp)) - 8.*exp(-(ak.*D).*(2.*trise+tp)) + 8.*exp(-(ak.*D).*(3.*trise+tp)) - 8.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 8.*exp(-(ak.*D).*(3.*trise+3.*tp)) + 8.*exp(-(ak.*D).*(4.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(5.*trise+2.*tp)) + 8.*exp(-(ak.*D).*(5.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(5.*trise+4.*tp)) - 8.*exp(-(ak.*D).*(6.*trise+3.*tp)) - 8.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 4.*exp(-(ak.*D).*(7.*trise+4.*tp)) + 12.*(ak.*D).^3.*trise.^3 + 16.*(ak.*D).^3.*trise.^2.*tp + 4.*exp(-(ak.*D).*(Delta-trise-tp)) - 2.*exp(-(ak.*D).*(Delta-trise-2.*tp)) + 4.*exp(-(ak.*D).*(Delta+trise+tp)) + 4.*exp(-(ak.*D).*(Delta-2.*trise-tp)) - 2.*exp(-(ak.*D).*(Delta+trise+2.*tp)) + 4.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 4.*exp(-(ak.*D).*(Delta-3.*trise-tp)) - 4.*exp(-(ak.*D).*(Delta+3.*trise+tp)) + 4.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) + 4.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 4.*exp(-(ak.*D).*(Delta+3.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(Delta-3.*trise-3.*tp)) - 4.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 4.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) - 2.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) - 2.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) - 4.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) - 4.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) - 2.*exp(-(ak.*D).*(Delta+5.*trise+4.*tp)) - 2.*exp(-(ak.*D).*(Delta-5.*trise-4.*tp)) + 4.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 4.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) + 4.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) + 4.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 2.*exp(-(ak.*D).*(Delta+7.*trise+4.*tp)) - 2.*exp(-(ak.*D).*(Delta-7.*trise-4.*tp)) + 16)./((ak.*D).^4.*trise.^2); 
        Etmp(Ind.tcos1,:,:) = Etmp_tcos1(Ind.tcos1,:,:) ; 
    end
    % tcos and n=2
    if any(Ind.tcos2)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tcos2 = 1/4*(24*exp(-(ak.*D).*tp) - 24*exp(-(ak.*D).*trise) - 48.*exp(-2.*(ak.*D).*trise) - 72.*exp(-(ak.*D).*Delta) - 120.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 24.*exp(-(ak.*D).*(Delta+2.*trise)) + 24.*exp(-(ak.*D).*(Delta-2.*trise)) - 12.*exp(-(ak.*D).*(Delta-tp)) - 12.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) + 36.*exp(-(ak.*D).*(trise+2.*tp)) - 24.*exp(-(ak.*D).*(2.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) - 72.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(3.*trise+3.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+3.*tp)) + 36.*exp(-(ak.*D).*(5.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(4.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(5.*trise+3.*tp)) - 24.*exp(-(ak.*D).*(6.*trise+3.*tp)) + 48.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(6.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(7.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(8.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(7.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(8.*trise+5.*tp)) + 24.*exp(-(ak.*D).*(9.*trise+5.*tp)) - 24.*exp(-(ak.*D).*(9.*trise+6.*tp)) - 24.*exp(-(ak.*D).*(9.*trise+7.*tp)) + 24.*exp(-(ak.*D).*(10.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(11.*trise+6.*tp)) + 24.*exp(-(ak.*D).*(11.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(11.*trise+8.*tp)) - 24.*exp(-(ak.*D).*(12.*trise+7.*tp)) - 24.*exp(-(ak.*D).*(12.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(13.*trise+8.*tp)) + 76.*(ak.*D).^3.*trise.^3 + 96.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) - 18.*exp(-(ak.*D).*(Delta-trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 12.*exp(-(ak.*D).*(Delta-2.*trise-tp)) - 18.*exp(-(ak.*D).*(Delta+trise+2.*tp)) + 12.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) + 36.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) + 36.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+3.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-3.*trise-3.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) - 18.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) - 18.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+4.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(Delta-4.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) + 12.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) - 24.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) - 24.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+6.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-6.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+7.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-7.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+8.*trise+4.*tp)) + 12.*exp(-(ak.*D).*(Delta-8.*trise-4.*tp)) - 6.*exp(-(ak.*D).*(Delta+7.*trise+6.*tp)) - 6.*exp(-(ak.*D).*(Delta-7.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+8.*trise+5.*tp)) + 12.*exp(-(ak.*D).*(Delta-8.*trise-5.*tp)) - 12.*exp(-(ak.*D).*(Delta+9.*trise+5.*tp)) - 12.*exp(-(ak.*D).*(Delta-9.*trise-5.*tp)) + 12.*exp(-(ak.*D).*(Delta+9.*trise+6.*tp)) + 12.*exp(-(ak.*D).*(Delta-9.*trise-6.*tp)) + 12.*exp(-(ak.*D).*(Delta+9.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-9.*trise-7.*tp)) - 12.*exp(-(ak.*D).*(Delta+10.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-10.*trise-7.*tp)) - 6.*exp(-(ak.*D).*(Delta+11.*trise+6.*tp)) - 6.*exp(-(ak.*D).*(Delta-11.*trise-6.*tp)) - 12.*exp(-(ak.*D).*(Delta+11.*trise+7.*tp)) - 12.*exp(-(ak.*D).*(Delta-11.*trise-7.*tp)) - 6.*exp(-(ak.*D).*(Delta+11.*trise+8.*tp)) - 6.*exp(-(ak.*D).*(Delta-11.*trise-8.*tp)) + 12.*exp(-(ak.*D).*(Delta+12.*trise+7.*tp)) + 12.*exp(-(ak.*D).*(Delta-12.*trise-7.*tp)) + 12.*exp(-(ak.*D).*(Delta+12.*trise+8.*tp)) + 12.*exp(-(ak.*D).*(Delta-12.*trise-8.*tp)) - 6.*exp(-(ak.*D).*(Delta+13.*trise+8.*tp)) - 6.*exp(-(ak.*D).*(Delta-13.*trise-8.*tp)) + 72)./(3.*(ak.*D).^4.*trise.^2) ; 
        Etmp(Ind.tcos2,:,:) = Etmp_tcos2(Ind.tcos2,:,:) ; 
    end
    % tcos and n=3
    if any(Ind.tcos3)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tcos3 = 1/4*((48.*exp(-2.*(ak.*D).*trise) - 72.*exp(-Delta.*(ak.*D)) - 60.*exp(-2.*(ak.*D).*tp) - 120.*(ak.*D).*tp - 30.*exp(-Delta.*(ak.*D)).*exp(-2.*(ak.*D).*trise) - 30.*exp((2.*trise-Delta).*(ak.*D)) + 36.*exp(-Delta.*(ak.*D)).*exp(-2.*(ak.*D).*tp) + 36.*exp((2.*tp-Delta).*(ak.*D)) - 96.*exp(-2.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 48.*exp(-2.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 36.*exp(-4.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 72.*exp(-4.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 36.*exp(-4.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 24.*exp(-6.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 48.*exp(-6.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 24.*exp(-6.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 12.*exp(-8.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 24.*exp(-8.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 12.*exp(-8.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) + 6.*exp(-10.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 12.*exp(-10.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) + 6.*exp(-10.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) - 12.*(ak.*D).^2.*tp.^2 + 40.*(ak.*D).^3.*tp.^3 + 12.*(ak.*D).^2.*tp.^2.*exp(-Delta.*(ak.*D)) + 120.*(ak.*D).^3.*trise.*tp.^2 + 24.*(ak.*D).*tp.*exp(-(ak.*D).*trise) + 60.*exp(-Delta.*(ak.*D)).*exp(-2.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 60.*exp((2.*trise+2.*tp-Delta).*(ak.*D)) - 30.*exp(-Delta.*(ak.*D)).*exp(-2.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 30.*exp((2.*trise+4.*tp-Delta).*(ak.*D)) + 24.*exp(-Delta.*(ak.*D)).*exp(-4.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 24.*exp((4.*trise+2.*tp-Delta).*(ak.*D)) - 48.*exp(-Delta.*(ak.*D)).*exp(-4.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 48.*exp((4.*trise+4.*tp-Delta).*(ak.*D)) + 24.*exp(-Delta.*(ak.*D)).*exp(-4.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 24.*exp((4.*trise+6.*tp-Delta).*(ak.*D)) - 18.*exp(-Delta.*(ak.*D)).*exp(-6.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 18.*exp((6.*trise+4.*tp-Delta).*(ak.*D)) + 36.*exp(-Delta.*(ak.*D)).*exp(-6.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 36.*exp((6.*trise+6.*tp-Delta).*(ak.*D)) - 18.*exp(-Delta.*(ak.*D)).*exp(-6.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 18.*exp((6.*trise+8.*tp-Delta).*(ak.*D)) + 12.*exp(-Delta.*(ak.*D)).*exp(-8.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 12.*exp((8.*trise+6.*tp-Delta).*(ak.*D)) - 24.*exp(-Delta.*(ak.*D)).*exp(-8.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 24.*exp((8.*trise+8.*tp-Delta).*(ak.*D)) + 12.*exp(-Delta.*(ak.*D)).*exp(-8.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) + 12.*exp((8.*trise+10.*tp-Delta).*(ak.*D)) - 6.*exp(-Delta.*(ak.*D)).*exp(-10.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 6.*exp((10.*trise+8.*tp-Delta).*(ak.*D)) + 12.*exp(-Delta.*(ak.*D)).*exp(-10.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) + 12.*exp((10.*trise+10.*tp-Delta).*(ak.*D)) - 6.*exp(-Delta.*(ak.*D)).*exp(-10.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) - 6.*exp((10.*trise+12.*tp-Delta).*(ak.*D)) + 6.*(ak.*D).^2.*tp.^2.*exp(-8.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) + 6.*(ak.*D).^2.*tp.^2.*exp(-12.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((trise-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-(ak.*D).*trise) - 24.*(ak.*D).*tp.*exp(-(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) - 24.*(ak.*D).*tp.*exp(-3.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 24.*(ak.*D).*tp.*exp(-3.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) + 24.*(ak.*D).*tp.*exp(-5.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) - 24.*(ak.*D).*tp.*exp(-5.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) - 24.*(ak.*D).*tp.*exp(-7.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) + 24.*(ak.*D).*tp.*exp(-7.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp(-9.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp(-9.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp(-11.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp(-11.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) - 6.*(ak.*D).^2.*tp.^2.*exp(-Delta.*(ak.*D)).*exp(-12.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) - 6.*(ak.*D).^2.*tp.^2.*exp((12.*trise+12.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp((trise+2.*tp-Delta).*(ak.*D)) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-3.*(ak.*D).*trise).*exp(-2.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp((3.*trise+2.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-3.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((3.*trise+4.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-5.*(ak.*D).*trise).*exp(-4.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((5.*trise+4.*tp-Delta).*(ak.*D)) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-5.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp((5.*trise+6.*tp-Delta).*(ak.*D)) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-7.*(ak.*D).*trise).*exp(-6.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp((7.*trise+6.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-7.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((7.*trise+8.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-9.*(ak.*D).*trise).*exp(-8.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((9.*trise+8.*tp-Delta).*(ak.*D)) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-9.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp((9.*trise+10.*tp-Delta).*(ak.*D)) + 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-11.*(ak.*D).*trise).*exp(-10.*(ak.*D).*tp) - 12.*(ak.*D).*tp.*exp((11.*trise+10.*tp-Delta).*(ak.*D)) - 12.*(ak.*D).*tp.*exp(-Delta.*(ak.*D)).*exp(-11.*(ak.*D).*trise).*exp(-12.*(ak.*D).*tp) + 12.*(ak.*D).*tp.*exp((11.*trise+12.*tp-Delta).*(ak.*D)) + 60)./(3.*(ak.*D).^4.*tp.^2)) ; 
        Etmp(Ind.tcos3,:,:) = Etmp_tcos3(Ind.tcos3,:,:) ; 
    end
    % tsin and n=1
    if any(Ind.tsin1)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tsin1 = 1/4*(24.*exp(-(ak.*D).*tp) - 24.*exp(-(ak.*D).*trise) - 12.*exp(-2.*(ak.*D).*trise) - 36.*exp(-Delta.*(ak.*D)) - 48.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 6.*exp(-(ak.*D).*(Delta+2.*trise)) + 6.*exp(-(ak.*D).*(Delta-2.*trise)) - 12.*exp(-(ak.*D).*(Delta-tp)) - 12.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) - 24.*exp(-(ak.*D).*(2.*trise+tp)) - 12.*exp(-(ak.*D).*(2.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(4.*trise+2.*tp)) + 16.*(ak.*D).^3.*trise.^3 + 24.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 12.*exp(-(ak.*D).*(Delta-2.*trise-tp)) + 12.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) + 6.*exp(-(ak.*D).*(Delta+2.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-2.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) + 6.*exp(-(ak.*D).*(Delta+4.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-4.*trise-2.*tp)) + 36) ./ (3.*(ak.*D).^4.*trise.^2) ; 
        Etmp(Ind.tsin1,:,:) = Etmp_tsin1(Ind.tsin1,:,:) ; 
    end
    % tsin and n=2
    if any(Ind.tsin2)
        trise = pulse.trise ; tp = pulse.tp ; trise = repmat(trise,[1 NR, Nk]) ; tp = repmat(tp,[1 NR, Nk]) ; 
        Etmp_tsin2 = 1/4*(48.*exp(-(ak.*D).*tp) - 24.*exp(-(ak.*D).*trise) - 36.*exp(-2.*(ak.*D).*trise) - 60.*exp(-Delta.*(ak.*D)) - 96.*(ak.*D).*trise + 12.*exp(-(ak.*D).*(Delta-trise)) + 12.*exp(-(ak.*D).*(Delta+trise)) + 18.*exp(-(ak.*D).*(Delta+2.*trise)) + 18.*exp(-(ak.*D).*(Delta-2.*trise)) - 24.*exp(-(ak.*D).*(Delta-tp)) - 24.*exp(-(ak.*D).*(Delta+tp)) - 24.*exp(-(ak.*D).*(trise+tp)) - 72.*exp(-(ak.*D).*(2.*trise+tp)) - 36.*exp(-(ak.*D).*(2.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(3.*trise+tp)) + 24.*exp(-(ak.*D).*(3.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+tp)) + 48.*exp(-(ak.*D).*(4.*trise+2.*tp)) + 24.*exp(-(ak.*D).*(4.*trise+3.*tp)) - 24.*exp(-(ak.*D).*(5.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(5.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(6.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(6.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(6.*trise+4.*tp)) + 24.*exp(-(ak.*D).*(7.*trise+3.*tp)) + 24.*exp(-(ak.*D).*(7.*trise+4.*tp)) - 12.*exp(-(ak.*D).*(8.*trise+4.*tp)) + 32.*(ak.*D).^3.*trise.^3 + 48.*(ak.*D).^3.*trise.^2.*tp + 12.*exp(-(ak.*D).*(Delta-trise-tp)) + 12.*exp(-(ak.*D).*(Delta+trise+tp)) + 36.*exp(-(ak.*D).*(Delta-2.*trise-tp)) + 36.*exp(-(ak.*D).*(Delta+2.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-tp)) + 18.*exp(-(ak.*D).*(Delta+2.*trise+2.*tp)) + 18.*exp(-(ak.*D).*(Delta-2.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-tp)) - 12.*exp(-(ak.*D).*(Delta+3.*trise+2.*tp)) - 12.*exp(-(ak.*D).*(Delta-3.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+tp)) - 24.*exp(-(ak.*D).*(Delta+4.*trise+2.*tp)) - 24.*exp(-(ak.*D).*(Delta-4.*trise-2.*tp)) - 12.*exp(-(ak.*D).*(Delta+4.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-4.*trise-3.*tp)) + 12.*exp(-(ak.*D).*(Delta+5.*trise+2.*tp)) + 12.*exp(-(ak.*D).*(Delta-5.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+5.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-5.*trise-3.*tp)) + 6.*exp(-(ak.*D).*(Delta+6.*trise+2.*tp)) + 6.*exp(-(ak.*D).*(Delta-6.*trise-2.*tp)) + 12.*exp(-(ak.*D).*(Delta+6.*trise+3.*tp)) + 12.*exp(-(ak.*D).*(Delta-6.*trise-3.*tp)) + 6.*exp(-(ak.*D).*(Delta+6.*trise+4.*tp)) + 6.*exp(-(ak.*D).*(Delta-6.*trise-4.*tp)) - 12.*exp(-(ak.*D).*(Delta+7.*trise+3.*tp)) - 12.*exp(-(ak.*D).*(Delta-7.*trise-3.*tp)) - 12.*exp(-(ak.*D).*(Delta+7.*trise+4.*tp)) - 12.*exp(-(ak.*D).*(Delta-7.*trise-4.*tp)) + 6.*exp(-(ak.*D).*(Delta+8.*trise+4.*tp)) + 6.*exp(-(ak.*D).*(Delta-8.*trise-4.*tp)) + 60) ./ (3.*(ak.*D).^4.*trise.^2) ;  
        Etmp(Ind.tsin2,:,:) = Etmp_tsin2(Ind.tsin2,:,:) ; 
    end
        
    %% final signal attenuation
    G = repmat(pulse.G, [1 NR]) ; 
    E = exp( - 2 .* (pulse.gamma*G).^2 .* sum(Bk .* Etmp,3)) ; 

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Bk, ak] = load_Bk_ak(R, structure) 
%% get Bk and ak 
    NR = length(R) ; 
    
    switch structure.geometry
        case 'plane'
            % R is the half distance between two planes
            krange = 1:30 ;             Nk = length(krange) ; 
            krange = reshape(krange,[1 1 Nk]) ; krange = repmat(krange, [1 NR 1]) ; 
            R = repmat(R, [1 1 Nk]) ; 
            Bk = 8.*(2*R).^2./pi^4./(2.*krange-1).^4 ; 
            ak = pi^2.*(2.*krange-1).^2./(2*R).^2 ; 

        case 'cylinder'
            muk = [ 1.8412, 5.3315, 8.5364, 11.7061, 14.8636, 18.0156, 21.1644, 24.3114, 27.4571, 30.6020, 33.7462, 36.8900, 40.0335, 43.1767, 46.3196, 49.4624, 52.6051, 55.7476, 58.8901, 62.0324, 65.1747, 68.3169, 71.4590, 74.6011, 77.7432, 80.8852, 84.0272, 87.1692, 90.3111  ]' ; 
            Nk = length(muk) ; 
            muk = reshape(muk, [1 1 Nk]) ; muk = repmat(muk,[1 NR 1]) ;             R = repmat(R,[1 1 Nk]) ; 
           Bk = 2.*(R./muk).^2./(muk.^2-1) ; 
           ak = (muk./R).^2 ; 

        case 'sphere'
            muk = [ 2.082,  5.941, 9.206, 12.405, 15.580, 18.743, 21.900, 25.053, 28.204, 31.353, 34.500, 37.646, 40.792, 43.937, 47.082, 50.226, 53.370, 56.514, 59.657, 62.801, 65.944, 69.087, 72.229, 75.372, 78.515, 81.657, 84.800, 87.942, 91.085  ]' ; 
            Nk = length(muk) ; 
            muk = reshape(muk, [1 1 Nk]) ; muk = repmat(muk,[1 NR 1]) ;             R = repmat(R,[1 1 Nk]) ; 
            Bk = 2.*(R./muk).^2./(muk.^2-2) ; 
            ak = (muk./R).^2 ; 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Bk, ak] = load_Bk_ak_hollow(Rin, Rout) 
    %% get Bk and ak for sphericalShell or hollowSphere
        lambdas = get_hollow_sphere_roots(Rin, Rout) ; 
        NumLambda = size(lambdas,1) ; 
        ak = lambdas.^2 ; 
%            Bkrange = 2.*Rin^3*Rout^3/(Rout^3-Rin^3)./lambdarange.^2 .* (sphbesselj_prim(1,lambdarange.*Rin)-sphbesselj_prim(1,lambdarange.*Rout)).^2 ...
%                ./(Rin^3.*(lambdarange.^2.*Rout^2-2).*sphbesselj_prim(1,lambdarange.*Rin).^2 - Rout^3.*(lambdarange.^2.*Rin^2-2).*sphbesselj_prim(1,lambdarange.*Rout).^2) ; 

%             Bkrange = 2./(Rout^3-Rin^3)./lambdarange.^2./Rout^3 .* (Rout^3.*sphbesselj1_core_prim(lambdarange.*Rin)-Rin^3.*sphbesselj1_core_prim(lambdarange.*Rout)).^2 ...
%                 ./ ((lambdarange.^2.*Rout^2-2).*sphbesselj1_core_prim(lambdarange.*Rin).^2 - (lambdarange.^2.*Rin^2-2).*sphbesselj1_core_prim(lambdarange.*Rout).^2) ; 
% 
        Rin = repmat(Rin, [NumLambda, 1]) ; 
        Rout = repmat(Rout, [NumLambda, 1]) ; 
        ratiotemp = Rout.^3./Rin.^3 .* sphbesselj1_core_prim(lambdas.*Rin) ./ sphbesselj1_core_prim(lambdas.*Rout) ; 
        Bk = 2*Rin.^3.*Rout.^3./(Rout.^3-Rin.^3)./lambdas.^2 .* (ratiotemp - 1).^2 ...
                ./ (Rin.^3.*(lambdas.^2.*Rout.^2-2).*ratiotemp.^2 - Rout.^3.*(lambdas.^2.*Rin.^2-2)) ; 

        % permute dimension
        ak = permute(ak, [3 2 1]) ; 
        Bk = permute(Bk, [3 2 1]) ; 
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lambdan = get_hollow_sphere_roots(a,b)
    %% get lambda for hollow sphere structure
    xstep = 0.00001 ;     xmax = 100 ; 
    x = xstep:xstep:xmax ; x = x' ; 
    numRoots = 20 ; 
    lambdan = zeros([numRoots, size(a,2)]) ; 
    
    y1 = sphbessely1_core_prim(x.*a).*sphbesselj1_core_prim(x.*b) - sphbesselj1_core_prim(x.*a).*sphbessely1_core_prim(x.*b) ; 
    % y1 = sphbessely_prim(1,x.*a).*sphbesselj_prim(1,x.*b) - sphbesselj_prim(1,x.*a).*sphbessely_prim(1,x.*b) ; 
    x = x + xstep ; 
    y2 = sphbessely1_core_prim(x.*a).*sphbesselj1_core_prim(x.*b) - sphbesselj1_core_prim(x.*a).*sphbessely1_core_prim(x.*b) ; 
    %y2 = sphbessely_prim(1,x.*a).*sphbesselj_prim(1,x.*b) - sphbesselj_prim(1,x.*a).*sphbessely_prim(1,x.*b) ; 

    for n = 1:size(a,2)
        root_tmp = x(y1(:,n).*y2(:,n)<0) ; 
        lambdan(:,n) = root_tmp(1:numRoots) ; 
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spherical bessel functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = sphbesselj1_core_prim(x)
    y = -2.*sin(x) + 2.*x.*cos(x) + x.^2.*sin(x) ; 
end

function y = sphbessely1_core_prim(x)
    y = 2.*cos(x) + 2.*x.*sin(x) - x.^2.*cos(x) ; 
end

function y = sphbesselj(n,x)
    y = besselj(n+0.5,x)./sqrt(x) ; 
end

function y = sphbessely(n,x)
    y = bessely(n+0.5,x)./sqrt(x) ; 
end

function y = sphbesselj_prim(n,x)
    %y = (n.*sphbesselj(n-1,x) - (n+1).*sphbesselj(n+1,x)) ./ (2*n+1) ; 
    y = n./x.*sphbesselj(n,x) - sphbesselj(n+1,x) ; 
end

function y = sphbessely_prim(n,x)
    %y = (n.*sphbessely(n-1,x) - (n+1).*sphbessely(n+1,x)) ./ (2*n+1) ; 
    y = n./x.*sphbessely(n,x) - sphbessely(n+1,x) ; 
end




