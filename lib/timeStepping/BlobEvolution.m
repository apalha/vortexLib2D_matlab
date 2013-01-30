function newBlobPositions = BlobEvolution(oldBlobs,deltaT,sigma,varargin)
% BlobEvolution advance vortex blobs in time.
%
%   newBlobPositions = BlobEvolution(oldBlobs,varargin)
%
%   Where:
%
% INPUTS
%
%   oldBlobs :: positions and circulations of the vortex blobs at the start
%               of the time step. oldBlobs(1,k) contains the x coordinate
%               of blob k. oldBlobs(2,k) contains the y coordinate
%               of blob k. oldBlobs(3,k) contains the circulation of blob k.
%               (type: real, dimension: [3 nBlobs])
%   deltaT :: the time step size
%             (type: real, dimension; [1])
%   sigma :: the square of radius of the blobs
%             (type: real, dimension: [1])
%
%   VARARGS
%       'TimeScheme' :: specify the integration scheme:
%                       'rk4' :: Runge-Kutta 4th order (default)
%                       'euler' :: explicit Euler 1st order
%                       (type: string, dimension: the size of the string)
%       'vInf' :: specify free stream velocity, default 0
%                 (type: real, dimension: [1 2])
%       'blockSize' :: specify the block size use in GPU computation,
%                      default 256
%                      (type: int, dimension: [1])
%
% OUTPUTS
%   newBlobPositions :: positions of the vortex blobs at the end of the
%                       time step. newBlobPositions(1,k) contains the x
%                       coordinate of blob k. newBlobPositions(1,k) contains
%                       the y coordinate of blob k.
%                       (type: real, dimension [2 nBlobs])
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/04 $

    % check for optional inputs
    
    % TimeScheme
    if any(strcmp('TimeScheme',varargin))
        optIndex = find(strcmp('TimeScheme',varargin));
        % if a time scheme is specified use it
        timeScheme = varargin{optIndex+1};
        if strcmp(timeScheme,'rk4')
            timeScheme = 0;
        end
        if strcmp(timeScheme,'euler')
            timeScheme = 1;
        end
    else
        % if no time scheme is specified, use rk4, Runge-Kutta 4th order
        timeScheme = 0;
    end
    
    % free stream velocity
    if any(strcmp('vInf',varargin))
        optIndex = find(strcmp('vInf',varargin));
        % if a free stream velocity is specified use it
        vInf = varargin{optIndex+1};
        vInfFlag = true;
    else
        % if no free stream velocity is specified, set it to zero
        vInfFlag = false;
    end
    
    % GPU block size
    if any(strcmp('blockSize',varargin))
        optIndex = find(strcmp('blockSize',varargin));
        % if a block size is specified use it
        blockSize = varargin{optIndex+1};
    else
        % if no block size is specified use 256
        blockSize = 256;
    end
    
    % check the number of particles and tranform it into a number of
    % particles multiple of the blocksize
    
    nBlobs = size(oldBlobs,2); % the number of particles
    
    nGhostParticles = ceil(nBlobs/blockSize)*blockSize-nBlobs; % the number of ghost particles
    
    oldBlobs = [oldBlobs zeros([3, nGhostParticles])]; % generate the new particles with extra ghost particles
    
    if timeScheme == 0
        % Runge-Kutta 4th order
        oldBlobsTempStage = oldBlobs;
        % compute the first stage
        k1 = (deltaT/6.0)*vv2parCC_mex_066DPGauss(oldBlobs,sigma,1.0/(2.0*pi),blockSize,eps);
        % compute the second stage
        oldBlobsTempStage(1:2,:) = oldBlobs(1:2,:) + k1*3; 
        k2 = (deltaT/3.0)*vv2parCC_mex_066DPGauss(oldBlobsTempStage,sigma,1.0/(2.0*pi),blockSize,eps);
        % compute the third stage
        oldBlobsTempStage(1:2,:) = oldBlobs(1:2,:) + k2*(3.0/2.0); 
        k3 = (deltaT/3.0)*vv2parCC_mex_066DPGauss(oldBlobsTempStage,sigma,1.0/(2.0*pi),blockSize,eps);
        % compute the fourth stage
        oldBlobsTempStage(1:2,:) = oldBlobs(1:2,:) + k3*3.0; 
        k4 = (deltaT/6.0)*vv2parCC_mex_066DPGauss(oldBlobsTempStage,sigma,1.0/(2.0*pi),blockSize,eps);
        
        % put the four stages together
        if vInfFlag
            % with free stream velocity
            newBlobPositions = oldBlobs(1:2,:) + k1 + k2 + k3 + k4 + deltaT*repmat(vInf,[1 size(oldBlobs,2)]);
        else
            % without free stream velocity
            newBlobPositions = oldBlobs(1:2,:) + k1 + k2 + k3 + k4;
        end
    elseif timeScheme == 1
        % Explicit Euler
        newBlobPositions = oldBlobs(1:2,:) + deltaT*vv2parCC_mex_066DPGauss(oldBlobs,sigma,1.0/(2.0*pi),blockSize,eps) + deltaT*repmat(vInf,[1 size(oldBlobs,2)]);
    end
    
    % return only the valid blobs
    newBlobPositions = newBlobPositions(1:2,1:nBlobs);
end