%% Add path to auxiliary functions
addpath('../../lib/inducedVelocities/')
addpath('../../lib/inducedVorticities/')
addpath('../../lib/initializationParticles/')
addpath('../../lib/plotting/')
addpath('../../lib/populationControl/')
addpath('../../lib/redistributionParticles/')
addpath('../../lib/velocityKernels/')
addpath('../../lib/vorticityKernels/')
addpath('../../lib/timeStepping/')

%% Definition of initial vorticity distribution

wExactFunction = @(x,y) wBlob(x,y,[0 1],50,0.5*sqrt(pi/12))... % Gaussian vortex located at [-1 0], with circulation (Gamma) equal to 50 and width (sigma) equal to 0.5*sqrt(pi/12)
                        + wBlob(x,y,[0 -1],-50,0.5*sqrt(pi/12)); % Gaussian vortex located at [1 0], with circulation (Gamma) equal to 50 and width (sigma) equal to 0.5*sqrt(pi/12)

%% Definition of the initial domain and blob parameters

blobsParameters.xBounds = [-10 10];     % x coordinates domain
blobsParameters.yBounds = [-10 10];     % y coordinates domain

blobsParameters.nBlobs = [128 128];     % number of blobs in both x and y directions

blobsParameters.overlap = 0.5;          % the overlap ratio between blobs

blobsParameters.populationControl = [1e6 1e6]*eps; % defines the minimum value
                                                    % of circulation to consider
                                                    % values below this are
                                                    % discarded and the
                                                    % maximul total
                                                    % circulation change
                                                    % when controlling
                                                    % population

blobsParameters.vInfinity = [0.0; 0.0];             % the velocity at infinity
                                                    
%% Definition of time stepping parameters

timeStepParametersBlobs.deltaT = 0.01;          % the step size in seconds
timeStepParametersBlobs.ntSteps = 2000;          % the number of time steps
timeStepParametersBlobs.stepDisplay = 10;      % the number of time steps between plot of vorticity
timeStepParametersBlobs.stepRedistribute = 10;  % the number of time steps between blob redistribution
timeStepParametersBlobs.stepPopulationControl = 10;  % the number of time steps between blob population control
timeStepParametersBlobs.timeScheme = 'rk4';     % time step scheme:
                                                %       'rk4'   :: Runge-Kutta 4th order
                                                %       'euler' :: basic forward Euler (testing purposes only)

                                                
%% Definition of GPU related parameters

gpuParameters.blockSize = 256;  % the size of the memory block in the GPU

%% Definition of the plotting domain

plotParameters.xBounds = [-5 30];  % x coordinates of plot domain
plotParameters.yBounds = [-5 5];  % y coordinates of plot domain

plotParameters.nPlotPoints = [128 128]; % number of plot point in both x and y directions

%% Compute initial blob properties and coordinates

% compute the spacing between the blobs
blobsParameters.deltaX = (blobsParameters.xBounds(2) - blobsParameters.xBounds(1))/blobsParameters.nBlobs(1);
blobsParameters.deltaY = (blobsParameters.yBounds(2) - blobsParameters.yBounds(1))/blobsParameters.nBlobs(2);

% determine if the spacing is equal in both directions, if not exit
if (blobsParameters.deltaX - blobsParameters.deltaY) >= 10*eps
    fprintf('\nBlob spacing not equal in x and y directions:\n    deltaX = %f\n    deltay = %f\n', blobsParameters.deltaX,blobsParameters.deltaY)
end

% compute blob properties
blobParameters.sigma = blobsParameters.deltaX/blobsParameters.overlap;   % the width of the blobs
blobParameters.area = blobsParameters.deltaX*blobsParameters.deltaY;     % the area of the cell used to blobify initial vorticity

% generate x and y coordinates of blobs
% recall that the blobs are cell centered (blobs close to boundary
% are at a distance 0.5*deltaXBlob of the boundary)
xBlob = (0:(blobsParameters.nBlobs(1)-1))'*blobsParameters.deltaX + 0.5*blobsParameters.deltaX + blobsParameters.xBounds(1);
yBlob = (0:(blobsParameters.nBlobs(2)-1))'*blobsParameters.deltaY + 0.5*blobsParameters.deltaY + blobsParameters.yBounds(1);

[blobsParameters.xGrid, blobsParameters.yGrid] = meshgrid(xBlob,yBlob);


%% Compute plot grid properties

% compute the plot grid spacings
plotParameters.deltaX = (plotParameters.xBounds(2)-plotParameters.xBounds(1))/(plotParameters.nPlotPoints(1)-1);
plotParameters.deltaY = (plotParameters.yBounds(2)-plotParameters.yBounds(1))/(plotParameters.nPlotPoints(2)-1);

% compute the x and y plot grid coordinates
xPlot = (0:(plotParameters.nPlotPoints(1)-1))'*plotParameters.deltaX + plotParameters.xBounds(1);
yPlot = (0:(plotParameters.nPlotPoints(2)-1))'*plotParameters.deltaY + plotParameters.yBounds(1);

[plotParameters.xBlobPlotGrid, plotParameters.yBlobPlotGrid] = meshgrid(xPlot,yPlot);


%% Compute exact vorticity

% compute the exact vorticity to plot
wExactPlot = wExactFunction(plotParameters.xBlobPlotGrid, plotParameters.yBlobPlotGrid);

% compute the exact vorticity at blob centers
wExact = wExactFunction(blobsParameters.xGrid, blobsParameters.yGrid);


%% Blobify the exact initial vorticity distribution

% compute initial circulations
Gamma = wExact*blobParameters.area;
GammaCorrected = Gamma;

% Beale's correction for initial circulation distribution

% A = BealeCorrectionMatrix(blobsParameters.xGrid,blobsParameters.yGrid,blobParameters.sigma,blobParameters.area);
% 
% % iterate
% for k=1:10
%     GammaCorrected = Gamma(:) + GammaCorrected(:) - A*GammaCorrected(:);
% end

% define the blobs
blobs = [blobsParameters.xGrid(:)'; blobsParameters.yGrid(:)'; GammaCorrected(:)'];


%% Clear GPU
GPU_mex_reset()


%% Compute blob evolution

% allocate memory space for total circulation
totalCirculation = zeros([timeStepParametersBlobs.ntSteps+1 1]);
totalCirculation(1) = sum(blobs(3,:));

% allocate memory space for number of blobs
nBlobsTime = zeros([timeStepParametersBlobs.ntSteps+1 1]);
nBlobsTime(1) = size(blobs,2);

% perform population control
blobs = PopulationControl(blobs,blobsParameters.populationControl(1),blobsParameters.populationControl(2));

% plot vorticity
PlotVorticity(blobs,...
              blobParameters.sigma*ones([size(blobs,2) 1]),...
              plotParameters.xBlobPlotGrid,...
              plotParameters.yBlobPlotGrid,...
              'FigNumber',2,...
              'Axis',[plotParameters.xBounds plotParameters.yBounds],...
              'AxisEqual',true);

colorbar
caxis([-50 50])



for tStep=1:timeStepParametersBlobs.ntSteps
    fprintf('\n ============================================ \n Time step :: %d \n',tStep);
    fprintf('    Number of blobs   :: %d \n',size(blobs,2));
    fprintf('    Total circulation :: %f \n',sum(blobs(3,:)));
    fprintf('    Circulation variation :: %e \n',(sum(blobs(3,:))-totalCirculation(1))/totalCirculation(1));
    
    % move particles
    blobs(1:2,:) = BlobEvolution(blobs,...
                                 timeStepParametersBlobs.deltaT,...
                                 blobParameters.sigma,...
                                 'blockSize',gpuParameters.blockSize,...
                                 'TimeScheme',timeStepParametersBlobs.timeScheme,...
                                 'vInf',blobsParameters.vInfinity);
    
    % redistribute blobs
    if mod(tStep,timeStepParametersBlobs.stepRedistribute)==0
        blobs = ParticleRedistribution2D(blobs,...
                                         blobsParameters.xBounds,...
                                         blobsParameters.yBounds,...
                                         blobsParameters.nBlobs(1),...
                                         blobsParameters.nBlobs(2));
                                     
        %disp(sprintf('Particles after redistribution = %d\n',size(blobs,2)))
    end
    
    % perform population control
    if mod(tStep,timeStepParametersBlobs.stepPopulationControl)==0
        blobs = PopulationControl(blobs,blobsParameters.populationControl(1),blobsParameters.populationControl(2));
        %disp(sprintf('Particles after population control = %d\n',size(blobs,2)))
    end
    
    % save total circulation of current time step
    totalCirculation(tStep+1) = sum(blobs(3,:));
    % save number of blobs of current time step
    nBlobsTime(tStep+1) = size(blobs,2);
    
    if mod(tStep,timeStepParametersBlobs.stepDisplay)==0
        tic
        PlotBlobs(blobs,blobParameters.sigma*ones(size(blobs,2),1),...
                  'FigNumber',2,...
                  'Axis',[-10,50,-10,10],...
                  'AxisEqual',true,...
                  'CirculationZFlag',true,...
                  'Alpha',0.5,...
                  'ShowEdges',false,...
                  'CirculationZ',true,...
                  'ClearFigure',true);
        toc
        
        % plot vorticity
        PlotVorticity(blobs,...
                      blobParameters.sigma*ones([size(blobs,2) 1]),...
                      plotParameters.xBlobPlotGrid,...
                      plotParameters.yBlobPlotGrid,...
                      'FigNumber',3,...
                      'Axis',[plotParameters.xBounds plotParameters.yBounds],...
                      'AxisEqual',true);

        colorbar
        caxis([-50 50])
        pause(0.1)
        %print('-f1','-djpeg99','-r100',sprintf('/home/gorkiana/results/twoVorticesBlobified/vorticity%03d_Alt.jpg',tDisplay))
        %print('-f2','-djpeg','-r100',sprintf('/home/gorkiana/results/twoVorticesBlobified/blobs%04d.jpg',tDisplay))
        %tDisplay = tDisplay + 1;
                
    end
    
    %blobs(1:2,:) = blobs(1:2,:) + repmat(deltaT*[0.05;0],[1 size(blobs,2)]);
end

%vWall(:,:,tStep+1) = vv2parCC_mex_066DPGaussTarget(blobs,[xWall(:)'; yWall(:)'],sqrt(blobRadiusSquared),1.0/(2.0*pi),256,eps);

% [vxWall(:,tStep+1) vyWall(:,tStep+1)] = ComputeVelocity(blobs,sqrt(blobRadiusSquared),xWall(:),yWall(:));

% clf(1)
% clf(2)
% PlotBlobs(blobs(:,1:nVortices),sqrt(blobRadiusSquared)*ones([size(blobs,2) 1]),'FigNumber',2,'Axis',[-10,30,-10,10],'AxisEqual',true);
% PlotVorticity(blobs(:,1:nVortices),sqrt(blobRadiusSquared)*ones([size(blobs,2) 1]),xPlotGrid,yPlotGrid,'FigNumber',1,'Axis',[-10,30,-10,10],'AxisEqual',true);
% pause(0.1)
% print('-f1','-djpeg','-r100',sprintf('/home/gorkiana/results/twoVorticesBlobified/vorticity%04d.jpg',tDisplay))
% print('-f2','-djpeg','-r100',sprintf('/home/gorkiana/results/twoVorticesBlobified/blobs%04d.jpg',tDisplay))
%         

GPU_mex_reset()

figure(4)
semilogy((1:timeStepParametersBlobs.ntSteps)*timeStepParametersBlobs.deltaT,abs(totalCirculation(2:end)-totalCirculation(1)))

figure(5)
loglog((0:timeStepParametersBlobs.ntSteps),nBlobsTime)


% for tStep=1:ntSteps
%    quiver(xWall(1:20:end),yWall(1:20:end),vxWall(1:20:end,tStep),vyWall(1:20:end,tStep),10)
%    axis equal
%    %axis([-10 10 -10 10])
%    pause(0.01)
% end