function newParticles = ParticleRedistribution(oldParticles,xBounds,nGridCells)
% ParticleInterpolation redistributes particles into a grid
%
%   newParticles = ParticleRedistribution(oldParticles,xBounds,nGridParticles)
%
%   Where:
%
% INPUTS
%
%   oldParticles :: the coordinates and circulation of the distorted
%                   particles
%                   (type: real, dimension: [2 nOldParticles])
%   xBounds :: the x bounds of the domain grid where new particles
%              are to be placed
%              (type: real, dimension: [1 2])
%
%   nGridCells :: the number of grid cells where new particles can be
%                     placed
%                     (type: integer, dimension: [1])
%
% OUTPUTS
%   newParticles :: the coordinates and circulation of the new particles
%                    (type: real, dimension: [2 nNewParticles])
%
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/10/05 $
    
    % compute grid spacing
    deltaX = (xBounds(2)-xBounds(1))/nGridCells;
    
    % determine the number of old particles
    nOldParticles = size(oldParticles,2);

    % compute the left side index of the grid particle of the kernel window of 
    % each particle
    xLIndices = (ceil((oldParticles(1,:)'-2*deltaX-xBounds(1)+0.5*deltaX)/deltaX)-1);

    % generate all the particles indices of the kernel of the particle
    xKernelIndices = [xLIndices xLIndices+1 xLIndices+2 xLIndices+3];

    % compute the coordinates of the grid particles of the kernel
    xKernelCoordinates = xKernelIndices*deltaX + 0.5*deltaX + xBounds(1);

    % compute the interpolating kernel for each grid particle
    kernelValues = M4Prime((xKernelCoordinates-repmat(oldParticles(1,:)',[1 4]))/deltaX).*repmat(oldParticles(2,:)',[1 4]);

    % compute the interpolation

    % find min and max indices
    minIndex = min(xKernelIndices(:));

    % add -minIndex+1 so that all indices are positive
    xKernelIndices = xKernelIndices - minIndex + 1;

    % compute the max values of kernel indices
    maxIndex = max(xKernelIndices(:));

    % allocate memory space for the index matrix
    indexMatrix = zeros([maxIndex 1]);

    % allocate memory space for the final interpolation and grid particles
    % coordinates vectors
    newParticles = zeros([2 nOldParticles*4]);

    % loop over interpolation kernel and update particles values
    nUniqueGridParticles = 0;

    for k=1:length(xKernelIndices(:))
        if indexMatrix(xKernelIndices(k))>0
            newParticles(2,indexMatrix(xKernelIndices(k))) = newParticles(2,indexMatrix(xKernelIndices(k))) + kernelValues(k);
        else
            nUniqueGridParticles = nUniqueGridParticles + 1;
            indexMatrix(xKernelIndices(k)) = nUniqueGridParticles;
            newParticles(1,nUniqueGridParticles) = xKernelCoordinates(k);
            newParticles(2,nUniqueGridParticles) = kernelValues(k);
        end
    end
    
    % only return the the particles with non-zero circulation
    newParticles = newParticles(:,1:nUniqueGridParticles);


end