function [vx vy] = ComputeVelocity(blobs,radius,xPlotGrid,yPlotGrid,varargin)
% ComputeVelocity computes the velocity of the blobs at xPlotGrid and
% yPlotGrid
%
%   [vx vy] = ComputeVelocity(blobs,radius,xPlotGrid,yPlotGrid,varargin)
%
%   Where:
%
% INPUTS
%
%   blobs :: the x and y coordinates of the blobs and the circulation G
%            (type: real, dimension: [3 nBlobs])
%   radius :: the radius of the blobs
%             (type: real, dimension: [1])
%   xPlotGrid :: the x coordinates where to plot the vorticity, generated
%                by meshgrid
%                (type: real, dimension: [nXGridPoints nYGridPoints])
%   yPlotGrid :: the y coordinates where to plot the vorticity, generated
%                by meshgrid
%                (type: real, dimension: [nXGridPoints nYGridPoints])
%   
%   VARARGS
%       'Method' :: specify which method to use in the computation of the
%                   induced velocities:
%                   0 is for direct calculation with the CPU in Matlab
%                     (default)
%                   1 is for direct calculation with the GPU (not used yet)
%                   2 is for fmm calculation with the CPU
%                   3 is for fmm calculation with the GPU
%                   (type: int, dimension [1])
%        'k' :: the k parameter that appears in the exponential. Typically
%               this value is 2 (the default value).
%               (type: real, dimension [1])
%               
%                   
%
% OUTPUTS
%   vx,vy :: the x and y velocities at the points xPlotGrid and yPlotGrid
%            (type: real, dimension [nXGridPoints nYGridPoints])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/14 $

    % check for optional inputs
    
    % kind of computation 
    if any(strcmp('Method',varargin))
        optIndex = find(strcmp('Method',varargin));
        % if a method is given, use that method to compute velocity
        method = varargin{optIndex+1};
    else
        % use default direct calculation
        method = 0;
    end
    
    % the k value 
    if any(strcmp('k',varargin))
        optIndex = find(strcmp('k',varargin));
        % if a value for k is given, use that value
        k = varargin{optIndex+1};
    else
        % use default value of k
        k = 2;
    end
    
    if method == 0 % CPU Matlab direct calculation
        % compute the number of blobs
        nBlobs = size(blobs,2);
        
        % compute the velocity
        vx = zeros(size(xPlotGrid)); % allocate memory space for velocity
        vy = zeros(size(xPlotGrid)); % allocate memory space for velocity
        
        % loop over all blobs, compute the vorticity and superimpose them
        for k=1:nBlobs
            vx = vx + vBlobX(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius,k);
            vy = vy + vBlobY(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius,k);
        end
        
    elseif method == 1 % GPU direct calculation (not yet implemented)
        disp('GPU direct calculation not available yet try another option!')
        
    elseif method == 2 % CPU fmm calculation
        targetsSize = size(xPlotGrid);
        targets = [xPlotGrid(:); yPlotGrid(:)];
        ksigmasqr = k*radius*radius;
        velocities = InducedVelocitiesFMM2dcpu(blobs,targets,ksigmasqr);
        vx = reshape(velocities(1,:),targetsSize);
        vy = reshape(velocities(2,:),targetsSize);
        
    elseif method == 3 % GPU fmm calculation
        targetsSize = size(xPlotGrid);
        targets = [xPlotGrid(:); yPlotGrid(:)];
        ksigmasqr = k*radius*radius;
        velocities = InducedVelocitiesFMM2dgpu(blobs,targets,ksigmasqr);
        vx = reshape(velocities(1,:),targetsSize);
        vy = reshape(velocities(2,:),targetsSize);
        
    end
        
end
