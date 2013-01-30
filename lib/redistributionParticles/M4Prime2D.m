function interpolationWeight = M4Prime2D(x,y)
% M4Prime2D computes the interpolating kernel of the M4 smooth interpolator
% for 2D.
%
%   interpolationWeight = M4Prime2D(x)
%
%   Where:
%
% INPUTS
%
%   x :: the x coordinates of the points where to compute the interpolation
%        weights
%        (type: real, dimension: (any))
%
%   y :: the y coordinates of the points where to compute the interpolation
%        weights
%        (type: real, dimension: dim(x))
%
% OUTPUTS
%   interpolationWeight :: the interpolation weight
%                          (type: real, dimension: dim(x))
%
% Implements:
%   interpolationWeight(x,y) = M4'(x)*M4'(y)
%
% See also M4PRIME.
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.0 $  $ Date: 2012/10/04 $

    % generate the interpolating kernel by tensor products
    xWeights = M4Prime(x);
    yWeights = M4Prime(y);
    
    % combine them
    interpolationWeight = xWeights.*yWeights;

end