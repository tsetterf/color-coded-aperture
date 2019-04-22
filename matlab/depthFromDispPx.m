function z = depthFromDispPx(d,a,b,c,zD,sPx)
%DEPTHFROMDISPPX Get the depth [m] from disparity [px]
%   Inputs:
%       d           Disparity [px]
%       a,b,c,zD    Calibration parameters
%       sPx         Camera pixel size [m]
%   Outputs:
%       z           Depth (z-coord) from the CCD [m]
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

    z = zD./(d*sPx/b-1+a)-c;
end
