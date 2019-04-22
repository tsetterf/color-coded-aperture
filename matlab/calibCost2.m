function cost = calibCost2(y,tPts_C,pGc,sPx)
%CALIBCOST2 Cost function in calibration optimization
%   Inputs:
%       y       Parameters [cGx, cGy, fEff, gEff]' of the projection model
%       tPts_C  3xN Matrix of pt source positions in CCD coords [m]
%       pGc     2xN Matrix of image circle centroid locs in CCD coords [px]
%       sPx     Pixel size [m]
%   Outputs:
%       cost    Squared Euclidean distance of point projection
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

    cGx = y(1); cGy = y(2); fEff = y(3); gEff = y(4);
    
    cost = 0;
    for i = 1:size(tPts_C,2)
        p = pinholeProject(tPts_C(:,i),fEff,cGx,cGy,sPx,gEff);
        cost = cost + norm(p - pGc(:,i))^2;
    end

end
