function [vx vy] = ComputeVelocity(blobs,radius,xPlotGrid,yPlotGrid)
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
% OUTPUTS
%   vx,vy :: the x and y velocities at the points xPlotGrid and yPlotGrid
%            (type: real, dimension [nXGridPoints nYGridPoints])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/14 $

    % check for optional inputs
    
    
    % compute the number of blobs
    nBlobs = size(blobs,2);
    
%     % Define the velocity field of a blob vortex
%     vBlobX = @(x,y,pVortex,gamma,epsilon) -(gamma./(2*pi*(eps+((x-pVortex(1)).^2) + ((y-pVortex(2)).^2)))).*...
%                                         (y-pVortex(2)).*(1-exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)));
% 
%     vBlobY = @(x,y,pVortex,gamma,epsilon) (gamma./(2*pi*(eps+((x-pVortex(1)).^2) + ((y-pVortex(2)).^2)))).*...
%                                         (x-pVortex(1)).*(1-exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)));

    % compute the velocity
    vx = zeros(size(xPlotGrid)); % allocate memory space for velocity
    vy = zeros(size(xPlotGrid)); % allocate memory space for velocity

    % loop over all blobs, compute the vorticity and superimpose them
    for k=1:nBlobs
        vx = vx + vBlobX(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius);
        vy = vy + vBlobY(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius);
    end  
    
end