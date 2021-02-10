function [Signal, SignalIn, SignalEx] = FcnSimulateSignal(sim, this)
%% simulate DWI signals using multiple parametric distributions

    %% simulate signals
    S0In = 0 ; 
    if sim.flag.DinDist == 'y'
        SignalIn = 0 ; 
        for nd=1:length(sim.ds)
            for nDin = 1:length(sim.Dins)
                SignalIn = SignalIn + sim.vd(nd)*sim.vDin(nDin)*(sim.ds(nd)^this.structure.Ndim)*mati.Physics.RestrictedDWISignal([sim.ds(nd), sim.Dins(nDin)]', this) ; 
                S0In = S0In + sim.vd(nd)*sim.vDin(nDin)*(sim.ds(nd)^this.structure.Ndim) ; 
            end
        end
    else
        SignalIn = 0 ; 
        for nd=1:length(sim.ds)
            SignalIn = SignalIn + sim.vd(nd)*(sim.ds(nd)^this.structure.Ndim)*mati.Physics.RestrictedDWISignal([sim.ds(nd), sim.Din]', this) ; 
            S0In = S0In + sim.vd(nd)*(sim.ds(nd)^this.structure.Ndim) ; 
        end
    end
    SignalIn = sim.vin .* SignalIn / S0In ; 

    
    %% extra cellular
    SignalEx = 0 ; 
    for nDex=1:length(sim.Dexs)
        for nbetaex = 1:length(sim.betaexs)
            SignalEx = SignalEx + sim.vDex(nDex)*sim.vbetaex(nbetaex)*exp(-this.pulse.b .* (sim.Dexs(nDex)+sim.betaexs(nbetaex).*this.pulse.f)) ; 
         end
    end
    SignalEx = sim.vex .* SignalEx + sim.vfree*exp(-this.pulse.b*this.structure.Dfree) ; 
    
    
    %% total signals
    Signal = SignalIn +  SignalEx ; 


end