function cost = calibCost(x,tPts_C,pRc,pGc,pBc,sPx,ind)
%CALIBCOST Cost function in calibration optimization
%   Inputs:
%       x       Parameters [ a, c, zD, b ]' of the disparity model
%       tPts_C  3xN Matrix of pt source positions in CCD coords [m]
%       pRc     2xN Matrix of image circle centroid locs in CCD coords [px]
%       pGc     2xN Matrix of image circle centroid locs in CCD coords [px]
%       pBc     2xN Matrix of image circle centroid locs in CCD coords [px]
%       sPx     Pixel size [m]
%       ind     Index of disparity BG=1, BR=2, RG=3.
%   Outputs:
%       cost    Cost magnitude
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

    % Try to model disparity = p/(depth+q) + r
    a = x(1); c = x(2); zD = x(3); b = x(4); 
    cost = 0;
    for i = 1:size(tPts_C,2)
        
        % Depth
        depth = abs(tPts_C(3,i));
        
        % Measured disparities
        if ind == 1             % dispBG
            dispMeas = sign(pGc(2,i)-pBc(2,i))*norm((pBc(:,i)-pGc(:,i))*sPx);
        elseif ind == 2         % dispBR
            dispMeas = sign(pRc(1,i)-pBc(1,i))*norm((pBc(:,i)-pRc(:,i))*sPx);
        else                    % dispRG
            dispMeas = sign(pGc(2,i)-pRc(2,i))*norm((pRc(:,i)-pGc(:,i))*sPx); 
        end
        
        % Modeled disparities
        dispMod = dispFromDepth(depth,a,b,c,zD);
        
        % Cost off difference between modeled and measured disparities
        cost = cost + 1e10 * norm(dispMod - dispMeas)^2;
        
    end

end
