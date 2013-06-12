function velocities = InducedVelocitiesFMM2dcpu(blobs,targets,ksigmasqr)
% InducedVelocitiesFMM2dcpu computes induced velocities at target points using
% a fast multipole method on a cpu.
%
%   velocities = InducedVelocitiesFMM2dcpu(blobs,targets,ksigmasqr)
%
%   Where:
%
% INPUTS
%
%   blobs :: the coordinates and circulation of the particles
%            (type: real, dimension: [3 nParticles])
%
%   targets :: the coordinates of the points where to compute induced the
%              velocities
%              (type: real, dimension: [2 nTargets])
%
%   ksigmasqr :: equal to k*sigma*sigma, where sigma is the core size of
%                the blobs and k is the spread parameter used.
%                (type: real, dimension: [1])
%
% OUTPUTS
%   velocities :: the induced velocities at the target points
%                 (type: real, dimension: [2 nTargets])
%
%
% This code is based on fmm2dgpu which was developed by Anders Goude. See 
% Goude, A., Engblom, S., Adaptive fast multipole methods on the GPU,
% The Journal of Scientific Computing, Volume 63, Issue 3, March 2013.
%
%   Copyright 2013 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2013/06/11 $
 
    % convert blob parameters into the parameters of the fmm solver
    alpha = 1.2564312086261696770;
    xopt = sqrt(alpha*ksigmasqr);
    cutoff = 5*xopt;
    
    % generate the vortices and targets in the format of the fmm solver
    blobsImag = blobs(1,:) + 1i*blobs(2,:);
    targetsImag = targets(1,:) + 1i*targets(2,:);
    
    % compute induced velocities
    velocitiesFMMcpuTemp = fmm2dcpu(blobsImag,blobs(3,:),targetsImag,'ndirect',35,'smooth','oseen','xopt',xopt,'cutoff',cutoff,'pot',1)/(2*pi);
    
    % convert velocities into the correct format
    velocities = zeros([2,length(velocitiesFMMcpuTemp)]);
    velocities(1,:) = -imag(velocitiesFMMcpuTemp);
    velocities(2,:) = -real(velocitiesFMMcpuTemp);
end