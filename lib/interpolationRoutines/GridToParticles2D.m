function newParticles = GridToParticles2D(oldParticles,wGrid,xBounds,yBounds,nGridCellsX,nGridCellsY)
% GridToParticles redistributes circulations from into a 2D grid
%
%   newParticles = ParticleRedistribution(oldParticles,xBounds,yBounds,nGridCellsX,nGridCellsY)
%
%   Where:
%
% INPUTS
%
%   oldParticles :: the coordinates of the distorted particles
%                   (type: real, dimension: [2 nOldParticles])
%
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
%   $ Revision: 1.0 $  $ Date: 2012/10/18 $
    
    % compute parameters
    h = (xBounds(2)-xBounds(1))/nGridCellsX;
    hy = (yBounds(2)-yBounds(1))/nGridCellsY;

    if (hy-h)>10*eps
        fprintf('\nSpacing of grid particles not equal: hx=%f and hy=%f\n\n',h,hy)
        return
    end
    
    % determine the number of old particles
    nOldParticles = size(oldParticles,2);

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
    kernelValues = M4Prime2D((xKernelCoordinates-repmat(oldParticles(1,:)',[1 16]))/h,(yKernelCoordinates-repmat(oldParticles(2,:)',[1 16]))/h);

    % compute the interpolation

    % allocate memory space for the final interpolation and grid particles
    % coordinates vectors
    newParticles = [oldParticles(1:2,:); zeros([1 nOldParticles])];
    
    xKernelIndices = xKernelIndices + 1;
    yKernelIndices = yKernelIndices + 1;
    
    % loop over interpolation kernel and update particles values
    for k=1:nOldParticles
        for m=1:16
            if (xKernelIndices(k,m)>0) && (yKernelIndices(k,m)>0) && (xKernelIndices(k,m)<=nGridCellsX) && (yKernelIndices(k,m)<=nGridCellsY)
                newParticles(3,k) = newParticles(3,k) + wGrid(xKernelIndices(k,m),yKernelIndices(k,m))*kernelValues(k,m);
            
            end
        end
    end

end