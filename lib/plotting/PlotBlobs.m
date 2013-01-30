function figNumber = PlotBlobs(blobs,radius,varargin)
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
%
%   VARARGS
%       'FigNumber' :: specify the figure number where to plot
%                      (type: integer, dimension: [1])
%       'Axis' :: specify the axis bounds
%                 (type: real, dimension: [1 4])
%       'AxisEqual' :: specify is equal scale for both axis is to be used
%                      (type: logic, dimension: [1])
%       'PlotThreshold' :: specify the threshold of circulation, below
%                          which blobs are not plotted
%                          (type: real, dimension: [1]) 
%       'CAxis' :: specify the colormap axis for giving colors to the blobs
%                          (type: real, dimension: [1 2])
%       'CirculationZ' :: flag that specifies if circulation is
%                                     used to set the Z coordinate of the blob
%                                     (type: logic, dimension: [1])
%       'Alpha' :: set the transparency alpha value of the blobs
%                  (type: real, dimension: [1])
%       'ShowEdges' :: Show the circle of the blobs as black or not
%                      (type: logic, dimension: [1])
%       'ClearFigure' :: Clear the figure or not
%                      (type: logic, dimension: [1])
%
% OUTPUTS
%   figNumber :: the number of the figure
%                (type: int, dimension [1])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.2 $  $ Date: 2013/01/29 $

% revision 1.1:
%   - Added 'PlotThreshold' option
%
% revision 1.2:
%   - Added 'CAxis' option
%
% revision 1.3:
%   - Blob plotting done with patches instead of rectangles
%   - Added option for plotting blobs in 3d with height being the
%     circulation
%   - Added option for transparency to the blobs
%   - Added option to show the edges as black
%   - Added option to clear figure

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
    
    % transparency
    if any(strcmp('ClearFigure',varargin))
        optIndex = find(strcmp('ClearFigure',varargin));
        ClearFigure = varargin{optIndex+1};
        
        if ClearFigure % clear figure if user says so
            clf(figNumber)
        end
    end
    
    % transparency
    if any(strcmp('Alpha',varargin))
        optIndex = find(strcmp('Alpha',varargin));
        blobsAlpha = varargin{optIndex+1};
    else
        % blobs are opaque
        blobsAlpha = 1;
    end
    
    % use circulation as z coordinate of blobs
    if any(strcmp('CirculationZ',varargin))
        optIndex = find(strcmp('CirculationZ',varargin));
        CirculationZFlag = varargin{optIndex+1};
    else
        % blobs are all at z=0
        CirculationZFlag = false;
    end
    
    % use circulation as z coordinate of blobs
    if any(strcmp('ShowEdges',varargin))
        optIndex = find(strcmp('ShowEdges',varargin));
        ShowEdgesFlag = varargin{optIndex+1};
    else
        % blobs are all at z=0
        ShowEdgesFlag = false;
    end
    
    if ShowEdgesFlag
        EdgeColor = 'k';
    else
        EdgeColor = 'None';
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
    
    % plot threshold
    if any(strcmp('PlotThreshold',varargin))
        optIndex = find(strcmp('PlotThreshold',varargin));
        PlotThreshold = varargin{optIndex+1};
    else
        % PlotThreshold is zero
        PlotThreshold = 0;
    end
    
    % caxis
    if any(strcmp('CAxis',varargin))
        optIndex = find(strcmp('CAxis',varargin));
        CAxis = varargin{optIndex+1};
        CAxisFlag = true;
    else
        % PlotThreshold is zero
        CAxisFlag = false;
    end
    
    % compute the number of blobs
    nBlobs = length(blobs(1,:));
   
% old version, slower
%
%     % generate the colors associated to the circulation value
%     cmap = 'jet';
%     cmap=colormap(cmap);                      % Set colormap
%     if CAxisFlag
%         yy=linspace(CAxis(1),CAxis(2),size(cmap,1));  % Generate range of color indices that map to cmap
%         cm = spline(yy,cmap',blobs(3,:));                  % Find interpolated colorvalue
%         cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
%         cm(cm<0)=0;
%     else
%         if (max(blobs(3,:))-min(blobs(3,:)))>100*eps
%             yy=linspace(min(blobs(3,:)),max(blobs(3,:)),size(cmap,1));  % Generate range of color indices that map to cmap
%             cm = spline(yy,cmap',blobs(3,:));                  % Find interpolated colorvalue
%             cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
%             cm(cm<0)=0;
%         else
%             % give the same color to all of the blobs (black)
%             cm = zeros([3 length(blobs(1,:))]);
%         end
%     end
%     % loop over the particles and plot the blobs
%     for k=1:nBlobs
%         if abs(blobs(3,k)) > PlotThreshold
%             rectangle('Position',[blobs(1,k)-radius(k) blobs(2,k)-radius(k) 2*radius(k) 2*radius(k) ],'Curvature',[1 1],'EdgeColor',cm(:,k),'FaceColor',cm(:,k))
%         end
%     end
%     
    % set the renderer to opengl in order to use transparency
    oldRenderer = get(gcf, 'Renderer');
    set(gcf, 'Renderer', 'opengl');
    
    % plot the blobs out of polygonal patches with 32 sides
    
    % generate the points of the approximated circles
    q = linspace(0,2*pi,32)';
    C = [cos(q) sin(q)];
    
    N = size(blobs,2);
    
    % plot all the patches
    if CirculationZFlag % use circulation as z component
        patch(bsxfun(@plus,bsxfun(@times,0.5*radius',repmat(C(:,1),[1 N])),blobs(1,:)), bsxfun(@plus,bsxfun(@times,0.5*radius',repmat(C(:,2),[1 N])),blobs(2,:)),repmat(blobs(3,:),[size(C,1) 1]),blobs(3,:),'edgecolor',EdgeColor);
    else % all blobs are at z=0
        patch(bsxfun(@plus,bsxfun(@times,0.5*radius',repmat(C(:,1),[1 N])),blobs(1,:)), bsxfun(@plus,bsxfun(@times,0.5*radius',repmat(C(:,2),[1 N])),blobs(2,:)),blobs(3,:),'edgecolor',EdgeColor);
    end
    
    % set the colorbar axis
    if CAxisFlag % use the one given by the user
        set(gca,'Clim',CAxis);
    else % use the minimum and maximum values of the circulation
        set(gca,'Clim',[min(blobs(3,:)) max(blobs(3,:))]);
    end
    
    % set the transparency of the blobs
    alpha(blobsAlpha)
    
    % show the colorbar
    colorbar
    %axis equal;

    
    % set equal axis
    if axisEqual
        axis equal
    end
    
    % set the axis
    if axisFlag
        axis(axisBounds)
    end
    
    set(0, 'DefaultFigureRendererMode', 'manual')
    set(0,'DefaultFigureRenderer',oldRenderer)
end