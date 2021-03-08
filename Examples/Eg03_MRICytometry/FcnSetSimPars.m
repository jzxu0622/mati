function sim = FcnSetSimPars(sim)
%% construct MRICytometry simulation parameters 
% INPUTS:
%       flag
% OUTPUTS: 
%       sim
% -------------------------------------------------------------------
% written by Junzhong Xu, Jan 29, 2020
% 
% -------------------------------------------------------------------

    %% Gaussian distributed d
%     sim.ds = (0.5:0.5:25)' ; sim.Dins = (0.2:0.2:3.2)' ; sim.Dexs = (0.2:0.2:3.2)' ; sim.betaexs = (0:0.5:10)' ; 
    dmin = 0.5 ; dd = 0.5 ; dmax = 25 ; sim.ds = (dmin:dd:dmax)' ; 
    Dinmin = 0.2 ; dDin = 0.2 ; Dinmax = 3.2 ; sim.Dins = (Dinmin:dDin:Dinmax)' ; 
    Dexmin = 0.2 ; dDex = 0.2 ; Dexmax = 3.2 ; sim.Dexs = (Dexmin:dDex:Dexmax)' ; 
    betaexmin = 0 ; dbetaex = 0.5 ; betaexmax = 10 ; sim.betaexs = (betaexmin:dbetaex:betaexmax)' ; 

    switch lower(sim.distribution)
        case 'gaussian'
            % ----------------- parameter distributions -----------------
            sim.vd= normpdf(sim.ds, sim.dcen, sim.dsigma) ; sim.vd= sim.vd/sum(sim.vd) ; sim.dmean = sum(sim.vd.* sim.ds)/sum(sim.vd) ; 
            sim.dvar = sum((sim.ds-sim.dmean).^2.*sim.vd)./sum(sim.vd) ; sim.dcov = sim.dsigma/sim.dmean ; 
            sim.vdvw = sim.vd.*sim.ds.^sim.Ndim ; sim.vdvw = sim.vdvw / sum(sim.vdvw) ; sim.dmeanvw = sum(sim.vd.* sim.ds.^(sim.Ndim))/sum(sim.vd.*sim.ds.^sim.Ndim) ; 
            sim.vDin = normpdf(sim.Dins, sim.Dincen, sim.Dinsigma) ; sim.vDin = sim.vDin/sum(sim.vDin) ; 
            sim.Dinmean = sum(sim.vDin .* sim.Dins) ; 

            sim.vDex = normpdf(sim.Dexs, sim.Dexcen, sim.Dexsigma) ; sim.vDex = sim.vDex/sum(sim.vDex) ; 
            sim.Dexmean = sum(sim.vDex .* sim.Dexs) ; 
            sim.vbetaex = normpdf(sim.betaexs, sim.betaexcen, sim.betaexsigma) ; sim.vbetaex = sim.vbetaex/sum(sim.vbetaex) ; 
            sim.betaexmean = sum(sim.vbetaex .* sim.betaexs) ; 

        case 'bimodal'
            % ----------------- parameter distributions -----------------
            sim.vd1 = normpdf(sim.ds, sim.dcen1, sim.dsigma1) ; 
            sim.vd2 = normpdf(sim.ds, sim.dcen2, sim.dsigma2) ; 
            sim.vd = sim.frac1*sim.vd1 + (1-sim.frac1)*sim.vd2 ; sim.vd= sim.vd/sum(sim.vd) ; sim.dmean = sum(sim.vd.* sim.ds) ; 
            sim.vdvw = sim.vd.*sim.ds.^sim.Ndim ; sim.vdvw = sim.vdvw / sum(sim.vdvw) ; sim.dmeanvw = sum(sim.vd.* sim.ds.^(sim.Ndim+1))/sum(sim.vd.*sim.ds.^sim.Ndim) ; 
            sim.dvar = sum((sim.ds-sim.dmean).^2.*sim.vd)./sum(sim.vd) ; sim.dcov = sqrt(sim.dvar)/sim.dmean ; 
            sim.vDin = normpdf(sim.Dins, sim.Dincen, sim.Dinsigma) ; sim.vDin = sim.vDin/sum(sim.vDin) ; 
            sim.Dinmean = sum(sim.vDin .* sim.Dins) ; 

            sim.vDex = normpdf(sim.Dexs, sim.Dexcen, sim.Dexsigma) ; sim.vDex = sim.vDex/sum(sim.vDex) ; 
            sim.Dexmean = sum(sim.vDex .* sim.Dexs) ; 
            sim.vbetaex = normpdf(sim.betaexs, sim.betaexcen, sim.betaexsigma) ; sim.vbetaex = sim.vbetaex/sum(sim.vbetaex) ; 
            sim.betaexmean = sum(sim.vbetaex .* sim.betaexs) ; 

        case 'gamma'
            sim.vd = flipud(gampdf(sim.ds-4, sim.dalpha, sim.dbeta)) ; 
            sim.vd= sim.vd/sum(sim.vd) ; sim.dmean = sum(sim.vd.* sim.ds)/sum(sim.vd) ; 
            sim.dvar = sum((sim.ds-sim.dmean).^2.*sim.vd)./sum(sim.vd) ; sim.dcov = sqrt(sim.dvar)/sim.dmean ; 
            sim.vdvw = sim.vd.*sim.ds.^sim.Ndim ; sim.vdvw = sim.vdvw / sum(sim.vdvw) ; sim.dmeanvw = sum(sim.vd.* sim.ds.^(sim.Ndim+1))/sum(sim.vd.*sim.ds.^sim.Ndim) ; 
            sim.vDin = normpdf(sim.Dins, sim.Dincen, sim.Dinsigma) ; sim.vDin = sim.vDin/sum(sim.vDin) ; 
            sim.Dinmean = sum(sim.vDin .* sim.Dins) ; 

            sim.vDex = normpdf(sim.Dexs, sim.Dexcen, sim.Dexsigma) ; sim.vDex = sim.vDex/sum(sim.vDex) ; 
            sim.Dexmean = sum(sim.vDex .* sim.Dexs) ; 
            sim.vbetaex = normpdf(sim.betaexs, sim.betaexcen, sim.betaexsigma) ; sim.vbetaex = sim.vbetaex/sum(sim.vbetaex) ; 
            sim.betaexmean = sum(sim.vbetaex .* sim.betaexs) ; 

        otherwise 
            error('the d distribution should be either gaussian or bimodal') ; 
    end
    
    %% adding free water to the simulation spectra
    if sim.vex~=0
        sim.vDex(end) = sim.vDex(end) + sim.vfree/sim.vex ; 
        sim.vbetaex(1) = sim.vbetaex(1) + sim.vfree/sim.vex ; 
    end
    
end

