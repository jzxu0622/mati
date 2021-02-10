function [ fitopts ] = CreateDictionary(this)
%% create Dictionary for MRI cytometry fitting
% 
% structure.modelName =  
%         'dNone_DinDist_DexNone'                %Dictionary is [1:NDin], i.e., diffusion distribution
%         'dDist_DinFixed_DexNone'                % Dictionary is [1:Nd]
%         'dDist_DinDist_DexNone'                   % Dictionary is [NDin * Nd]
%         'dDist_DinFixed_DexDist'                   % Dictionary is [1:Nd, 1:Dex]
%         'dDist_DinDist_DexDist'                      % Dictionary is [NDin * Nd, 1:Dex]
%         'dDist_DinFixed_DexDisper'               % Dictionary is [1 * Nd, Nbetaex * NDex]            
%         'dDist_DinDist_DexDisper'                  % Dictionary is [NDin * Nd, Nbetaex * NDex]
%
% -------------------------------------------------------------------------------------------------------------------------

    %% preliminary
    structure = this.structure ; pulse = this.pulse ; fitopts = this.defaultFitopts ; 
    
    %% fitting parameters
    switch structure.modelName
        case 'dNone_DinDist_DexNone'               %Dictionary is [1:NDin]
            fitopts.ds = [] ; fitopts.Nd = 0 ; 
            fitopts.Dins= (fitopts.Dinmin:fitopts.DDin:fitopts.Dinmax)' ;  fitopts.NDin = length(fitopts.Dins) ; 
            fitopts.Dexs = [] ; fitopts.NDex = 0 ; 
            fitopts.betaexs = [] ; fitopts.Nbetaex = 0 ; 
        case 'dDist_DinFixed_DexNone'                 % Dictionary is [1:Nd]
            fitopts.ds= (fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;  fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins = structure.Din ; fitopts.NDin = 1 ; 
            fitopts.Dexs = [] ; fitopts.NDex = 0 ; 
            fitopts.betaexs = [] ; fitopts.Nbetaex = 0 ; 
        case 'dDist_DinDist_DexNone'                      % Dictionary is [NDin * Nd]
            fitopts.ds=(fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;  fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins= (fitopts.Dinmin:fitopts.DDin:fitopts.Dinmax)' ; fitopts.NDin = length(fitopts.Dins) ; 
            fitopts.Dexs = [] ; fitopts.NDex = 0 ; 
            fitopts.betaexs = [] ; fitopts.Nbetaex = 0 ; 
        case 'dDist_DinFixed_DexDist'                     % Dictionary is [1:Nd, 1:Dex]
            fitopts.ds = (fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;              fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins = structure.Din ; fitopts.NDin = 1 ; 
            fitopts.Dexs= (fitopts.Dexmin:fitopts.DDex:fitopts.Dexmax)' ; fitopts.NDex = length(fitopts.Dexs) ; 
            fitopts.betaexs = 0 ; fitopts.Nbetaex = 1 ; 
        case 'dDist_DinDist_DexDist'                  % Dictionary is [NDin * Nd, 1:Dex]
            fitopts.ds=(fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;  fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins= (fitopts.Dinmin:fitopts.DDin:fitopts.Dinmax)' ; fitopts.NDin = length(fitopts.Dins) ; 
            fitopts.Dexs= (fitopts.Dexmin:fitopts.DDex:fitopts.Dexmax)' ; fitopts.NDex = length(fitopts.Dexs) ; 
            fitopts.betaexs = 0 ; fitopts.Nbetaex = 1 ; 
        case 'dDist_DinFixed_DexDisper'                        % Dictionary is [1 * Nd, Nbetaex * NDex]            
            fitopts.ds = (fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;  fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins = structure.Din ; fitopts.NDin = 1 ; 
            fitopts.Dexs = (fitopts.Dexmin:fitopts.DDex:fitopts.Dexmax)' ; fitopts.NDex = length(fitopts.Dexs) ; 
            fitopts.betaexs= (fitopts.betaexmin:fitopts.Dbetaex:fitopts.betaexmax)' ; fitopts.Nbetaex = length(fitopts.betaexs) ; 
        case 'dDist_DinDist_DexDisper'                        % Dictionary is [NDin * Nd, Nbetaex * NDex]
            fitopts.ds = (fitopts.dmin:fitopts.Dd:fitopts.dmax)' ;  fitopts.Nd = length(fitopts.ds) ; 
            fitopts.Dins = (fitopts.Dinmin:fitopts.DDin:fitopts.Dinmax)' ; fitopts.NDin = length(fitopts.Dins) ; 
            fitopts.Dexs = (fitopts.Dexmin:fitopts.DDex:fitopts.Dexmax)' ; fitopts.NDex = length(fitopts.Dexs) ; 
            fitopts.betaexs= (fitopts.betaexmin:fitopts.Dbetaex:fitopts.betaexmax)' ; fitopts.Nbetaex = length(fitopts.betaexs) ; 
        otherwise
            error('%s: The input structure.modelName is not recognized',mfilename) ; 
    end
    % -------------------------- indices -----------------------------
    indAll = false([fitopts.NDin*fitopts.Nd+fitopts.Nbetaex*fitopts.NDex,1]) ; 
    fitopts.indIn = indAll ; fitopts.indIn(1:(fitopts.NDin*fitopts.Nd)) = true ; fitopts.Nin = sum(fitopts.indIn(:)) ; 
    fitopts.indEx = indAll ; fitopts.indEx((fitopts.NDin*fitopts.Nd+1):end) = true ; fitopts.Nex = sum(fitopts.indEx(:)) ; 
    % ---------------------- create Dictionary -------------------
    fitopts.Dictionary = zeros([pulse.Nacq, length(indAll)]) ; 
    [X_d, Y_Din] = meshgrid(fitopts.ds, fitopts.Dins) ; 
    [X_Dex, Y_betaex] = meshgrid(fitopts.Dexs, fitopts.betaexs) ; 
    
    % cell-volumn-weighted Dictionary
    for n = 1:fitopts.Nin
        fitopts.Dictionary(:,n) = mati.Physics.RestrictedDWISignal([X_d(n), Y_Din(n)]', this) ;  
    end
    for n = 1:fitopts.Nex
        fitopts.Dictionary(:,fitopts.Nin+n) = exp(-pulse.b .* (X_Dex(n)+Y_betaex(n).*pulse.f)) ; % Dex = Dex0 + beta*f 
    end
    
    % non-cell-volumn-weighted Dictionary2
    for n = 1:fitopts.Nin
        fitopts.Dictionary2(:,n) = X_d(n)^structure.Ndim*mati.Physics.RestrictedDWISignal([X_d(n), Y_Din(n)]', this) ;  
    end
    fitopts.Dictionary2 = fitopts.Dictionary2 ./ (sum(fitopts.ds(:).^structure.Ndim)) ; 
    
end
