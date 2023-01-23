function [structure, pulse] = SetCoordRange(structure, pulse)
%% function to set grid points coordinates range
% Definition of coordinates: X: rows, Y: columns, Z: slices
%
%
% ----------------------------------------------------------------------------------------------------

    %% calculate angles
    % angle between gradient and z-axis
    pulse.gTheta = acos(dot([0 0 1], pulse.gdir)) ; 
    % angle between x-axis and the projection of gradient onto the x-y plane
    pulse.gPhi = acos(dot([1 0 0], [pulse.gdir(1), pulse.gdir(2), 0])) ; 
    
    %% coordinate range along gradient direction
    switch structure.Ndim
        case 1  % 1D
            if abs(pulse.gTheta-pi/2)>eps || abs(pulse.gPhi)>eps
                pulse.gTheta = pi/2 ; 
                pulse.gPhi = 0 ; 
                warning('DWI:GdirectionWrong', '1D tissue, gradient should be along x-axis! \n set pulse.gTheta=pi/2, pulse.gPhi=0.') ; 
            end
            structure.coordRange = ((1:structure.N(1))' - (structure.N(1)+1)/2) .* structure.dr(1) ; 
        case 2  % 2D
            if abs(pulse.gTheta-pi/2)>eps 
                pulse.gTheta = pi/2 ; 
                warning('DWI:GdirectionWrong', '2D tissue, gradient should be on x-y plane! \n set pulse.gTheta=pi/2.') ; 
            end
            [X,Y] = meshgrid(1:structure.N(1), 1:structure.N(2)) ; 
            X = (X - (structure.N(1)+1)/2) .* structure.dr(1) ; 
            Y = (Y - (structure.N(2)+1)/2) .* structure.dr(2) ; 
            structure.coordRange = X.*cos(pulse.gPhi) + Y.*sin(pulse.gPhi) ; 
            % swith x and y axis, because matlab set horizontal as x, vertical as y
            structure.coordRange = structure.coordRange' ; 
            structure.coordRange = reshape(structure.coordRange, structure.size, 1) ; 
        case 3  % 3D
            [X,Y,Z] = meshgrid(1:structure.N(1), 1:structure.N(2), 1:structure.N(3)) ; 
            X = (X - (structure.N(1)+1)/2) .* structure.dr(1) ; 
            Y = (Y - (structure.N(2)+1)/2) .* structure.dr(2) ; 
            Z = (Z - (structure.N(3)+1)/2) .* structure.dr(3) ; 
            structure.coordRange = X.*sin(pulse.gTheta).*cos(pulse.gPhi) + Y.*sin(pulse.gTheta).*sin(pulse.gPhi) + Z.*cos(pulse.gTheta) ; 
            % swith x and y axis, because matlab set horizontal as x, vertical as y
            coordRangeTmp = zeros(structure.N) ; 
            for i = 1:structure.N(3)
                coordRangeTmp(:,:,i) = structure.coordRange(:,:,i)' ; 
            end
            structure.coordRange = reshape(coordRangeTmp, structure.size, 1) ; 
            clear coord_range_tmp
    end

    %% convenient factors for speeding up 
    % structure.phase_factor = 1i .* pulse.gamma .* pulse.dt .* structure.coord_range ;  
    structure.rpbcFactorX = 1i .* pulse.gamma .* structure.N(1) .* structure.dr(1) .* sin(pulse.gTheta) .* cos(pulse.gPhi) ; 
    structure.rpbcFactorY = 1i .* pulse.gamma .* structure.N(2) .* structure.dr(2) .* sin(pulse.gTheta) .* sin(pulse.gPhi) ; 
    structure.rpbcFactorZ = 1i .* pulse.gamma .* structure.N(3) .* structure.dr(3) .* cos(pulse.gTheta) ; 


end


