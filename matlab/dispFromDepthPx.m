function d = dispFromDepthPx(depth,a,b,c,zD,sPx)
%DISPFROMDEPTHPX Get the disparity [px] from depth [m]
%   Inputs:
%       depth       Depth (z-coord) from the CCD [m]
%       a,b,c,zD    Calibration parameters
%       sPx         Camera pixel size [m]
%   Outputs:
%       d           Disparity at the input depth [px]
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

    d = dispFromDepth(depth,a,b,c,zD) / sPx;
    
end
