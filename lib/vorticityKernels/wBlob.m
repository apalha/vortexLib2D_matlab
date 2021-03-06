function w = wBlob(x,y,pVortex,gamma,epsilon)
% wBlob computes the vorticity distribution generated by vortex
% of circulation gamma, located at pVortex and core spreading epsilon

    w = gamma*(exp(-(((x-pVortex(1)).^2) + ((y-pVortex(2)).^2))/(2*epsilon*epsilon)))/(2*pi*epsilon*epsilon);
end
