function d = dispFromDepth(depth,a,b,c,zD)
%DISPFROMDEPTH Get the disparity [m] from depth [m]
%   Inputs:
%       depth       Depth (z-coord) from the CCD [m]
%       a,b,c,zD    Calibration parameters
%   Outputs:
%       d           Disparity at the input depth [m]
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

    d = b * ( 1 - a + zD./(depth + c) );
    
end
