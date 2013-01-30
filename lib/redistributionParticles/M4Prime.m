function interpolationWeight = M4Prime(x)
% M4Prime computes the interpolating kernel of the M4 smooth interpolator.
%
%   interpolationWeight = M4Prime(x)
%
%   Where:
%
% INPUTS
%
%   x :: the points where to compute the interpolation weights
%        (type: real, dimension: any)
%
% OUTPUTS
%   interpolationWeight :: the interpolation weight
%                          (type: real, dimension: dim(x))
%
% Implements:
%             /
%            | 0                           if |x| > 2
%            |
%  M4'(x) =<  0.5*((2-|x|)^2)*(1-|x|)     if |x| \in [1,2]
%            |
%            | 1 - 2.5*(x^2) + 1.5*(|x|^3) if |x| \in [0,1]
%             \
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/09/27 $

    % allocate memory space for interpolationWeight
    interpolationWeight = zeros(size(x));
    
    % compute each condition
    
    % |x| \in [1,2]
    selection = (abs(x)>1) & (abs(x)<2); % compute the values in this range
    interpolationWeight(selection) = 0.5*((2-abs(x(selection))).^2).*(1-abs(x(selection))); % evaluate the interpolation weights
    
    % |x| \in [0,1]
    selection = (abs(x)>=0) & (abs(x)<1); % compute the values in this range
    interpolationWeight(selection) = 1 - 2.5*(x(selection).^2) + 1.5*(abs(x(selection)).^3); % evaluate the interpolation weights
end