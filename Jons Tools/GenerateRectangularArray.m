function [arraymat] = GenerateRectangularArray(m,n,deltax,deltaz)
    % Creates a 3D matrix of position vectors for an mxn planar array
    % on the y = 0 plane

    x = deltax*0.5*ones(m,1)*linspace(-1,1,n);
    z = deltaz*0.5*(ones(n,1)*linspace(-1,1,m))';
    y = zeros(m,n);
    arraymat(:,:,1) = x';
    arraymat(:,:,2) = y';
    arraymat(:,:,3) = z';
end