function newParticles = ParticleRedistribution2D(oldParticles,xBounds,yBounds,nGridCellsX,nGridCellsY)
% ParticleInterpolation2D redistributes particles into a 2D grid
%
%   newParticles = ParticleRedistribution(oldParticles,xBounds,yBounds,nGridCellsX,nGridCellsY)
%
%   Where:
%
% INPUTS
%
%   oldParticles :: the coordinates and circulation of the distorted
%                   particles
%                   (type: real, dimension: [3 nOldParticles])
%   xBounds,yBounds :: the x and y bounds of the domain grid where new particles
%                      are to be placed
%                      (type: real, dimension: [1 2])
%
%   nGridCellsX,nGridcellsY :: the number of grid cells where new particles can be
%                              placed in x and y directions
%                     (type: integer, dimension: [1])
%
% OUTPUTS
%   newParticles :: the coordinates and circulation of the new particles
%                   (type: real, dimension: [3 nNewParticles])
%
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/10/05 $
    
    if isempty(oldParticles)
        newParticles = oldParticles;
    else
        % compute parameters
        h = (xBounds(2)-xBounds(1))/nGridCellsX;
        hy = (yBounds(2)-yBounds(1))/nGridCellsY;

        if (hy-h)>10*eps
            fprintf('\nSpacing of grid particles not equal: hx=%f and hy=%f\n\n',h,hy)
            return
        end

        % determine the number of old particles
        %nOldParticles = size(oldParticles,2);

        % compute the left side index and bottom side index of the grid particle of the kernel window of 
        % each particle
        xLIndices = (ceil((oldParticles(1,:)'-2*h-xBounds(1)+0.5*h)/h)-1);
        yBIndices = (ceil((oldParticles(2,:)'-2*h-yBounds(1)+0.5*h)/h)-1);

        % generate all the particles indices of the kernel of the particle
        xKernelIndices = [repmat(xLIndices,[1 4]) repmat(xLIndices+1,[1 4]) repmat(xLIndices+2,[1 4]) repmat(xLIndices+3,[1 4])];
        yKernelIndices = [yBIndices yBIndices+1 yBIndices+2 yBIndices+3 yBIndices yBIndices+1 yBIndices+2 yBIndices+3 yBIndices yBIndices+1 yBIndices+2 yBIndices+3 yBIndices yBIndices+1 yBIndices+2 yBIndices+3];

        % compute the coordinates of the grid particles of the kernel
        xKernelCoordinates = xKernelIndices*h + 0.5*h + xBounds(1);
        yKernelCoordinates = yKernelIndices*h + 0.5*h + yBounds(1);

        % compute the interpolating kernel for each grid particle
        kernelValues = M4Prime2D((xKernelCoordinates-repmat(oldParticles(1,:)',[1 16]))/h,(yKernelCoordinates-repmat(oldParticles(2,:)',[1 16]))/h).*repmat(oldParticles(3,:)',[1 16]);

        % compute the interpolation

        % find min and max indices
        minIndexX = min(xKernelIndices(:));
        minIndexY = min(yKernelIndices(:));

        % add -minIndexX+1 and -minIndexY+1 so that all indices are positive
        xKernelIndices = xKernelIndices - minIndexX + 1;
        yKernelIndices = yKernelIndices - minIndexY + 1;

        % compute the max values of kernel indices
        maxIndexX = max(xKernelIndices(:));
        maxIndexY = max(yKernelIndices(:));

        [newParticlesIndexI newParticlesIndexJ newParticlesData] = find(sparse(xKernelIndices(:),yKernelIndices(:),kernelValues(:),maxIndexX,maxIndexY));

        % only return the the particles with non-zero circulation
        newParticles = [(newParticlesIndexI'+minIndexX-1)*h+0.5*h+xBounds(1); (newParticlesIndexJ'+minIndexY-1)*h+0.5*h+yBounds(1); newParticlesData'];
    end
end