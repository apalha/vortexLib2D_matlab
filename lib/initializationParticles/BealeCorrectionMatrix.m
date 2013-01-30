 function A = BealeCorrectionMatrix(xBlobGrid,yBlobGrid,blobRadius,blobCellArea,varargin)
% BealeCorrectionMatrix generates Beale's correction matrix to improve blob
% initialization data.
%
%   Asparse = BealeCorrectionMatrix(xBlobGrid,yBlobGrid,blobRadius)
%
%   Where:
%
% INPUTS
%
%   xBlobGrid :: x coordinates of the center of the blobs
%                (type: real, dimension: [nBlobs 1])
%   yBlobGrid :: y coordinates of the center of the blobs
%                (type: real, dimension: [nBlobs 1])
%   blobRadius :: the radius of the blobs
%                 (type: real, dimension: [1])
% 
% Optional inputs:
%   
%   'wBlob' :: a function defining the blob vorticity kernel, if no
%              vorticity kernel is specified, then the default one is used
%
% OUTPUTS
%   A :: Beale's correction matrix
%        (type: real, dimension [nBlobs nBlobs] (sparse))
%
%   Copyright 2012 Artur Palha
%   $ Revision: 1.1 $  $ Date: 2013/01/28 $

%   $ Revision: 1.0 $  $ Date: 2012/09/03 $
%   $ Revision: 1.1 $  $ Date: 2013/01/28 $
%       - added optional wBlob fuction as input

    % vorticity kernel
    if any(strcmp('wBlob',varargin))
        optIndex = find(strcmp('wBlob',varargin));
        % if a vorticity kernel is given, use it
        wBlob = varargin{optIndex+1};
    end

    % compute the number of blobs
    nBlobs = length(xBlobGrid(:));
    
    % allocate memory space for A matrix
    A = zeros(nBlobs);
    
    for k=1:nBlobs
        A(:,k) = blobCellArea*reshape(wBlob(xBlobGrid,yBlobGrid,[xBlobGrid(k) yBlobGrid(k)],1,blobRadius),[nBlobs 1]);
    end
    
    % convert to a sparse matrix: all values smaller than 10*eps are
    % descarded
    nonZeroData = (A>(10e8)*eps); % find the nonzero values
    
    Adata = A(nonZeroData); % take only the nonzero data
    
    clear('A'); % get rid of the whole old matrix, which is not needed any more
    
    [jIndices,iIndices] = meshgrid(1:nBlobs); % generate all the indices of the data in the matrix
    
    jIndices = jIndices(nonZeroData); % take only the ones associated to the non-zero values
    iIndices = iIndices(nonZeroData); % take only the ones associated to the non-zero values
    
    
    A = sparse(iIndices,jIndices,Adata,nBlobs,nBlobs); % generate the sparse matrix
end
