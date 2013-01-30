function w = ComputeVorticity(blobs,radius,xPlotGrid,yPlotGrid)
% ComputeVorticity computes the vorticity of the blobs at xPlotGrid and
% yPlotGrid
%
%   w = ComputeVorticity(blobs,radius,xPlotGrid,yPlotGrid,varargin)
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
%   w :: the vorticity at the points xPlotGrid and yPlotGrid
%        (type: real, dimension [nXGridPoints nYGridPoints])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/13 $

    % check for optional inputs
    
    
    % compute the number of blobs
    nBlobs = size(blobs,2);
    
%     % define the vorticity function
%     wBlob = @(x,y,pVortex,gamma,epsilon) gamma*(exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)))/(2*pi*epsilon*epsilon);
    
    % compute the vorticity
    w = zeros(size(xPlotGrid)); % allocate memory space for vorticity
    
    % loop over all blobs, compute the vorticity and superimpose them
    for k=1:nBlobs
        w = w + wBlob(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius);
    end

    
end