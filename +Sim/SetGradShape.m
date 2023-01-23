function pulse = SetGradShape(pulse)
%% 


%% gradients series
% if structure.flag.useT2 == 'y'
%     pulse.gShape = zeros(1,round((pulse.TE)/pulse.dt)) ; 
% else
%     pulse.gShape = zeros(1,round((pulse.Delta+pulse.delta)/pulse.dt)) ; 
% end

pulse.gShape = zeros(1,round((pulse.TE)/pulse.dt)) ; 

switch lower(pulse.shape)
    case 'pgse'
        ntdelta = round(pulse.delta/pulse.dt) ; 
        ntDelta = round(pulse.Delta/pulse.dt) ; 
        pulse.gShape(1:ntdelta) = 1 ; 
        pulse.gShape((ntDelta+1):(ntDelta+ntdelta)) = -1 ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntdelta
%                 pulse.gShape(nt) = 1 ; 
%             elseif nt>ntDelta && nt<=(ntDelta + ntdelta)
%                 pulse.gShape(nt) = - 1 ; 
%             end
%         end

    case 'tpgse'
        ntp = round(pulse.tp/pulse.dt) ; 
        ntrise = round(pulse.trise/pulse.dt) ; 
%         nDelta = round(pulse.Delta/pulse.dt) ; 
        % 1st gradient
        g11 = (1:ntrise)/ntrise ; g12 = ones([1 ntp]) ; g13 = 1-(1:ntrise)/ntrise ; 
        g1 =[g11, g12, g13]  ; 
        pulse.gShape = g1 ; 
        % between two gradients
        if pulse.Delta - pulse.delta > pulse.dt 
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end
        % 2nd gradient
        pulse.gShape = [pulse.gShape, -g1] ; 
        % after 2nd gradiet and before TE
        if pulse.TE - pulse.Delta - pulse.delta > pulse.dt
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.TE - pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end
        
    case 'sin'
        ntdelta = round(pulse.delta/pulse.dt) ; 
        ntDelta = round(pulse.Delta/pulse.dt) ; 
        for nt = 1:length(pulse.gShape)
            if nt<=ntdelta
                pulse.gShape(nt) = sin(2*pi/pulse.T*(nt-1)*pulse.dt) ; 
            elseif nt>ntDelta && nt<=(ntDelta + ntdelta)
                pulse.gShape(nt) = - sin(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
            end
        end
  
%     case 'sin2'
%         % sin2, OGSE stops simulation right after the 1st gradient cycle
%         pulse.gShape = zeros(1,round((pulse.T)/pulse.dt)) ; 
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntsigma
%                 pulse.gShape(nt) = sin(2*pi/pulse.T*(nt-1)*pulse.dt) * square(pi/pulse.T*(nt-1)*pulse.dt); 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  sin(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)) * square(pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
%             end
%         end

    case 'cos'
        ntdelta = round(pulse.delta/pulse.dt) ; 
        ntDelta = round(pulse.Delta/pulse.dt) ; 
        for nt = 1:length(pulse.gShape)
            if nt<=ntdelta
                pulse.gShape(nt) = cos(2*pi/pulse.T*(nt-1)*pulse.dt) ; 
            elseif nt>ntDelta && nt<=(ntDelta + ntdelta)
                pulse.gShape(nt) = - cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
            end
        end
        
        
%     case 'cos2'
%         % cos2, OGSE stops simulation right after the 1st gradient cycle
%         pulse.gShape = zeros(1,round((pulse.T)/pulse.dt)) ; 
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntsigma
%                 pulse.gShape(nt) = pulse.g * cos(2*pi/pulse.T*(nt-1)*pulse.dt) * square(pi/pulse.T*2*(nt-1)*pulse.dt); 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  pulse.g * cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)) * square(pi/pulse.T*2*((nt-1)*pulse.dt-pulse.Delta)); 
%             end
%         end


    case 'tcos'
        ntp = round(pulse.tp/pulse.dt) ; 
        ntrise = round(pulse.trise/pulse.dt) ; 
        % gradient segments
        g1 = (1:ntrise)/ntrise ; g2 = ones([1 ntp]) ; g3 = 1-(1:ntrise)/ntrise ; g4 = -(1:ntrise)/ntrise ; g5 = -ones([1 ntp+floor(ntrise/2)]) ; g6 = -1+(1:ntrise)/ntrise ; g7 = (1:ntrise)/ntrise ; g8 = -g5 ; 
        g = [g1, g2, g3, g4,g5] ; 
        for nn=1:(pulse.n-1)
            g = [g, g5 g6 g7 g8 g8 g3 g4 g5] ; 
        end
        g = [g, g5, g6, g7, g2, g3] ; 
        % 1st gradient
        pulse.gShape = g ; 
        % between two gradients
        if pulse.Delta - pulse.delta > pulse.dt
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end
        % 2nd gradient
        pulse.gShape = [pulse.gShape, -g] ; 
        % after 2nd gradient & before TE
        if pulse.TE - pulse.Delta - pulse.delta > pulse.dt
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.TE - pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end        
        
    case 'tsin'
        ntp = round(pulse.tp/pulse.dt) ; 
        ntrise = round(pulse.trise/pulse.dt) ; 
        % gradient segments
        g1 = (1:ntrise)/ntrise ; g2 = ones([1 ntp]) ; g3 = 1-(1:ntrise)/ntrise ; g4 = -(1:ntrise)/ntrise ; g5 = -ones([1 ntp]) ; g6 = -1+(1:ntrise)/ntrise ; 
        g = [g1, g2, g3, g4,g5, g6] ; 
        for nn=1:(pulse.n-1)
            g = [g, g] ; 
        end
        % 1st gradient
        pulse.gShape = g ; 
        % between two gradients
        if pulse.Delta - pulse.delta > pulse.dt
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end
        % 2nd gradient
        pulse.gShape = [pulse.gShape, -g] ; 
        % after 2nd gradient & before TE
        if pulse.TE - pulse.Delta - pulse.delta > pulse.dt
            pulse.gShape = [pulse.gShape, zeros(1, round((pulse.Delta - pulse.delta)/pulse.dt))] ; 
        end
        
        
        %     case 'cos22'
%         % cos22, OGSE stops simulation right after the 1st gradient cycle
%         pulse.gShape = zeros(1,round((pulse.T)/pulse.dt)) ; 
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntsigma
%                 pulse.gShape(nt) = cos(2*pi/pulse.T*(nt-1)*pulse.dt) * square(pi/pulse.T*(nt-1)*pulse.dt); 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)) * square(pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
%             end
%         end
% 
% 
%     case 'square'
%         % square, OGSE stops simulation right after the 1st gradient cycle
%         pulse.gShape = zeros(1,round((pulse.T)/pulse.dt)) ; 
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntsigma
%                 pulse.gShape(nt) = square(2*pi/pulse.T*(nt-1)*pulse.dt) ; 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  square(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
%             end
%         end
%         
%         
%     case 'square2'
%         % square2
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=ntsigma
%                 pulse.gShape(nt) = square(2*pi/pulse.T*(nt-1)*pulse.dt) * square(pi/pulse.T*(nt-1)*pulse.dt); 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  square(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)) * square(pi/pulse.T*(nt-1)*pulse.dt); 
%             end
%         end
% 
% 
%     case 'sawtooth'
%         % sawtooth
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             nnn = floor((nt-1)*pulse.dt/pulse.T)*pulse.T+pulse.T/2 ; 
%             if nt<=ntsigma
%                 pulse.gShape(nt) = (1-tripuls((nt-1)*pulse.dt-nnn,pulse.T)*2) ; 
%             elseif nt>ntDelta && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) =  (1-tripuls((nt-1)*pulse.dt-nnn,pulse.T)*2); 
%             end
%         end

%     case 'apodised-cos'
%         % cos
%         ntsigma = round(pulse.delta/pulse.dt) ; 
%         ntDelta = round(pulse.Delta/pulse.dt) ; 
%         nT4 = round(pulse.T/4/pulse.dt) ; 
%         for nt = 1:length(pulse.gShape)
%             if nt<=nT4
%                 pulse.gShape(nt) = sign(cos(2*pi/pulse.T*(nt-1)*pulse.dt)) * abs(sin(4*pi/pulse.T*(nt-1)*pulse.dt)) ; 
%             elseif nt>nT4 && nt<=ntsigma-nT4
%                 pulse.gShape(nt) = cos(2*pi/pulse.T*(nt-1)*pulse.dt) ; 
%             elseif nt>ntsigma-nT4 && nt<=ntsigma
%                 pulse.gShape(nt) = sign(cos(2*pi/pulse.T*(nt-1)*pulse.dt)) * abs(sin(4*pi/pulse.T*(nt-1)*pulse.dt)) ; 
%             elseif nt>ntDelta && nt<=(ntDelta + nT4)
%                 pulse.gShape(nt) = - sign(cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta))) * abs(sin(4*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta))) ; 
%             elseif nt>ntDelta+nT4 && nt<=ntDelta+ntsigma-nT4
%                 pulse.gShape(nt) = - cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta)); 
%             elseif nt>(ntDelta + ntsigma - nT4) && nt<=(ntDelta + ntsigma)
%                 pulse.gShape(nt) = - sign(cos(2*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta))) * abs(sin(4*pi/pulse.T*((nt-1)*pulse.dt-pulse.Delta))) ; 
%             end
%         end

    otherwise
        error('pulse.shape should be "pgse", "tpgse", "sin", "cos", or "tcos"') ; 
end


%% integral of gradients
pulse.gShapeIntegral = cumsum(pulse.gShape) .* pulse.dt ; 

