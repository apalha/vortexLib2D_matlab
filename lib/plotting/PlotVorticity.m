function figNumber = PlotVorticity(blobs,radius,xPlotGrid,yPlotGrid,varargin)
% PlotBlobs plots vortex blobs with radius corresponding to the size of the
% blob and color corresponding to the circulation associated to the blob.
%
%   figNumber = PlotBlobs(blobs,h,varargin)
%
%   Where:
%
% INPUTS
%
%   blobs :: the x and y coordinates of the blobs and the circulation G
%            (type: real, dimension: [3 nBlobs])
%   radius :: the radius of the blobs
%             (type: real, dimension: [nBlobs 1])
%   xPlotGrid :: the x coordinates where to plot the vorticity, generated
%                by meshgrid
%                (type: real, dimension: [nXGridPoints nYGridPoints])
%   yPlotGrid :: the y coordinates where to plot the vorticity, generated
%                by meshgrid
%                (type: real, dimension: [nXGridPoints nYGridPoints])
%
%   VARARGS
%       'FigNumber' :: specify the figure number where to plot
%                      (type: integer, dimension: [1])
%       'Axis' :: specify the axis bounds
%                 (type: real, dimension: [1 4])
%       'AxisEqual' :: specify is equal scale for both axis is to be used
%                      (type: logic, dimension: [1])
%
% OUTPUTS
%   figNumber :: the number of the figure
%                (type: int, dimension [1])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/03 $

    % check for optional inputs
    
    % figure number
    if any(strcmp('FigNumber',varargin))
        optIndex = find(strcmp('FigNumber',varargin));
        figNumber = varargin{optIndex+1};
        % if a figure number is given, use that figure to plot
        figure(figNumber)
    else
        % generate a new figure
        figNumber = figure();
    end
    
    % axis values
    if any(strcmp('Axis',varargin))
        optIndex = find(strcmp('Axis',varargin));
        axisBounds = varargin{optIndex+1};
        % if axis are given, set the flag
        axisFlag = true;
    else
        % no axis specification is made
        axisFlag = false;
    end
    
    % axis equal
    if any(strcmp('AxisEqual',varargin))
        optIndex = find(strcmp('AxisEqual',varargin));
        axisEqual = varargin{optIndex+1};
    else
        % axis equal is not set
        axisEqual = false;
    end
    
    % number of contour levels
    if any(strcmp('nContourLevels',varargin))
        optIndex = find(strcmp('nContourLevels',varargin));
        nContourLevels = varargin{optIndex+1};
    else
        % if no number of contour levels is given use 100
        nContourLevels = 100;
    end
    
    % compute the number of blobs
    nBlobs = size(blobs,2);
    
    % define the vorticity function
    wBlob = @(x,y,pVortex,alpha,epsilon) alpha*(exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)))/(2*pi*epsilon*epsilon);
    
    % compute the corrected blobified vorticity
    wBlobified = zeros(size(xPlotGrid)); % allocate memory space for vorticity

    % loop over all vortices, compute the vorticity and superimpose them
    for k=1:nBlobs
        wBlobified = wBlobified + wBlob(xPlotGrid,yPlotGrid,[blobs(1,k) blobs(2,k)],blobs(3,k),radius(k));
    end
    max(wBlobified(:))
    % plot the vorticity
    %contourf(xPlotGrid,yPlotGrid,wBlobified,nContourLevels,'LineColor','none');
    %pcolor(xPlotGrid,yPlotGrid,wBlobified);
    surf(xPlotGrid,yPlotGrid,wBlobified);
    shading interp;
    
    % set equal axis
    if axisEqual
        axis equal
    end
    
    % set the axis
    if axisFlag
        axis(axisBounds)
    end
    
    
end