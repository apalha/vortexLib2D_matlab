function newParticles = PopulationControl(oldParticles,gThreshold,gRelativeThreshold)
% PopulationControl controls the number of particles by discarding
% particles with circulation higher than gThreshold
%
%   newParticles = PopulationControl(oldParticles,gThreshold,gRelativeThreshold)
%
%   Where:
%
% INPUTS
%
%   oldParticles :: the coordinates and circulation of the distorted
%                   particles
%                   (type: real, dimension: [3 nOldParticles])
%
%   gThreshold :: the minimum value of absolute circulation to be
%                 considered. Values lower than this are discarded.
%                 (type: real, dimension: [1])
%
%   gRelativeThreshold :: the maximum value of the ratio flaggedCirculation/totalCirculation
%                         that can be considered
%                         (type: real, dimension: [1])
%
% OUTPUTS
%   newParticles :: the coordinates and circulation of the new particles
%                   (type: real, dimension: [3 nNewParticles])
%
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/10/08 $
    
    % compute the total absolute circulation
    totalCirculation = sum(abs(oldParticles(3,:)));
    
    % determine the particles with circulation smaller then gThreshold
    % compute how much absolute circulation is discarded, if 
    % flaggedCirculation/totalCirculation < gRelativeThreshold
    % then discard particles, if not, decrease gThreshold by one order
    % of magnitude
    circulationConservationFlag = false; % initialize circulation flag that 
                                         % determines if the amount of particles
                                         % discarded is acceptable in terms
                                         % of circulation loss
    while ~circulationConservationFlag
        flaggedParticles = (abs(oldParticles(3,:))<gThreshold); % determine the particles to discard
        flaggedCirculation = sum(abs(oldParticles(3,flaggedParticles))); % determine the total flagged circulation
        if (flaggedCirculation/totalCirculation)<gRelativeThreshold
            circulationConservationFlag = true; % the circulation change is acceptable
        else
            gThreshold = gThreshold/10.0; % decrease the threshold
        end
    end
    
    % return the non-discarded particles
    newParticles = oldParticles(:,~flaggedParticles);
    


end