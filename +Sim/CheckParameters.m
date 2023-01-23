function [structure, pulse] = CheckParameters(structure, pulse)


%% ********************* test the consistancy of parameters ***************
% CFL condition, for explicit discretization scheme
if max(max(structure.D).*pulse.dt./structure.dr.^2) > 0.5/structure.Ndim
    error(['D*dt/dx^2=',num2str(structure.D*pulse.dt/max(structure.dr)^2),' too large!']) ; 
end

% q_dx effect
if max(pulse.gamma.*max(pulse.gShapeIntegral).*pulse.G.*structure.dr ./pi) > 0.101
    error(['q*dx/pi=', num2str(pulse.gamma * max(pulse.gShapeIntegral) *pulse.G* max(structure.dr)/pi),' too large!']) ; 
end

% OGSE gradient 
if strcmp(pulse.shape, 'sin') || strcmp(pulse.shape, 'cos') || strcmp(pulse.shape, 'tcos') || strcmp(pulse.shape, 'tsin')
    % OGSE period cannot be too short (less than 10*dt)
    if pulse.T < 10*pulse.dt
        disp(['T=', num2str(pulse.T/pulse.dt),' dt, too small']) ; 
    end
    
    if pulse.delta < pulse.T
        error('T should be smaller than sigma') ; 
    end
%     % OGSE: integer # of period inside each gradient
%     if abs(round(pulse.delta/pulse.T)*pulse.T - pulse.delta) > pulse.dt 
%         error('each gradient should have integer # of oscillations!') ; 
%     end
    
end


% check pulse.grad
if sum(pulse.gShape) > 1
    error('pulse.grad error! pulse sequence does not have two identical gradients on both side of pi pulse so the spins can not rephase completely. ') ; 
end




end

