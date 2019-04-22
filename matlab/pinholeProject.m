function p = pinholeProject(t,f,cX,cY,sPx,g)
%PINHOLEPROJECT Project a 3D point to a 2D image plane
%   Inputs:
%       t       3D point [x y z]' [px]
%       f       Focal length [m]
%       cX      Principal point x coordinate [px]
%       cY      Principal point y coordinate [px]
%       sPx     Pixel size [m]
%       g       Gap between front and rear principal planes (default 0) [m]
%   Outputs:
%       p       2D point [u v]' [px]
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Camera matrix for pinhole model
f = f/sPx;
K = [f 0 cX;
     0 f cY;
     0 0  1];

 if nargin < 6
    g = 0;
 end
 
% Rotate from conventional optical frame to conventional image frame
t = [-1 0 0; 0 1 0; 0 0 -1]*t;  

p = K * (t-[0 0 g]');
p = p(1:2)./p(3);

end

