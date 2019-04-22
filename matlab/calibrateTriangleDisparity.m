% This script loads the following files:
%   A calib.csv file, with format:
%       Range,PointId,ImagePSF,ImageAprilTags
%       1,511,DSC0324,DSC0325
%       ...
%   A camera.csv file, with format:
%       ImagePlanePointId,MountPointId,LensPointId,FocusReference,PixelSize
%       508,509,510,Inf,4.88e-6
%   A calib.idx file output from the total station. 
% The script loads the image _ImagePSF.JPG at each range and asks the user 
% to draw a box around the PSF. The images must contain a triangle 
% disparity PSF on a black background. All images should be taken at the 
% same reference focal distance zR. Determination of the filter crosstalk
% is performed separately in calibrateTriangleColorSpace.m. Previously 
% processed data can also be loaded from a .mat file to avoid the manual 
% steps.
% Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Image extension
ext = 'JPG';

% Camera name
camName = 'dalsa5100';

% Whether to load and save data
loadData = 1;
saveData = 0;

% Name of image at first depth to exclude from finding pos'ns pGcm and pBcm
imNameExcludeCentroids = 'DAL0013-short.JPG';

% Data folder
folder = '../data/calib-glass50-led-08.24.2018/';
% folder = '../data/calib-glass85-led-08.24.2018/';

% Focal length spec [m]
f = 50e-3;

% Whether calibration was completed
calibCompleted = 0;

% Not loading data -------------------------------------------------------
if loadData == 0

    % Get the calibration ids table
    calibTable = readtable([folder 'calib.csv'],'Delimiter',',');

    % Get the camera location table and pixel size
    camLocTable = readtable([folder 'camera.csv'],'Delimiter',',');
    zR          = camLocTable.FocusReference(1);
    sPx         = camLocTable.PixelSize(1);

    % Get the total station output
    readPts  = -1;
    tsFileId = fopen([folder 'calib.idx']);
    tsLine   = fgetl(tsFileId);     % total station file line
    tsPts    = [];                  % arr w/ rows [PointId East North Elev]
    while (ischar(tsLine))

        % Start reading points after string 'POINTS' and stop at 'THEMINFO'
        if contains(tsLine,'POINTS')
            readPts = 0;
        elseif contains(tsLine,'THEMINFO')
            break;
        end

        % Read points
        if readPts >= 0

            % Skip first point, but read other points
            if readPts >= 2
                tsLine = erase(tsLine,'"');
                lSplit = strsplit(tsLine,'\s*,\s*','DelimiterType', ...
                                        'RegularExpression');
                tsPts(end+1,:) = [str2double(lSplit{2})         ...
                    str2double(lSplit{3}) str2double(lSplit{4}) ...
                    str2double(lSplit{5})]';
            end

            readPts = readPts + 1;
        end

        % Get next line
        tsLine = fgetl(tsFileId);

    end
    fclose(tsFileId);

    % Camera point locations in the camera coordinate frame (origin at 
    % physical center of CCD, with +z from lens to CCD, +y down toward 
    % tripod mount) [m]
    % Nikkor Calib 2018.06.28 --------------------------------------------
    % tPo_C = [ [ -0.0290  -0.0585  0      ]' ... % p1, cam body target
    %           [ -0.00625  0.0640  0      ]' ... % p2, tripod mount target
    %           [ -0.0366   0      -0.0868 ]' ];  % p3, lens barrel target
    % Zeiss Glass f=50mm Calib 2018.08.24 --------------------------------
    tPo_C = [ [  0.02946   -0.00203  0        ]' ... % p1, cam body target
              [  0.01676    0.15990 -0.00592  ]' ... % p2, tripod mount tgt
              [  0.03683    0.00147 -0.07658  ]' ];  % p3, lens barrel tgt
    % Zeiss Glass f=85mm Calib 2018.08.24 --------------------------------
    % tPo_C = [ [ 0.02946  -0.00203  0        ]' ... % p1, cam body tgt
    %           [ 0.01676   0.15990 -0.00592  ]' ... % p2, tripod mount tgt
    %           [ 0.03980   0.00201 -0.08890  ]' ];  % p3, lens barrel tgt
    % ------------------------------------------------------------------------
    tPm_C = mean(tPo_C,2);                      % centroid location
    tP_C  = tPo_C - tPm_C;                      % refer all pts to centroid

    % Camera point locations in the total station coordinate frame
    nCamPoses = length(camLocTable.ImagePlanePointId);
    minAppPt    = zeros(nCamPoses,1); 
    maxAppPt    = zeros(nCamPoses,1);
    RTtoC       = zeros(3,3,nCamPoses);
    tCtoT_C     = zeros(3,nCamPoses);
    for i = 1:nCamPoses
        tPo_T = [ ...
          tsPts(tsPts(:,1)==camLocTable.ImagePlanePointId(i),2:4)'... % p1
          tsPts(tsPts(:,1)==camLocTable.MountPointId(i),2:4)'     ... % p2
          tsPts(tsPts(:,1)==camLocTable.FocusGripPointId(i),2:4)' ];  % p3
        tPm_T = mean(tPo_T,2);              % centroid location
        tP_T  = tPo_T - tPm_T;              % refer all points to centroid
        minAppPt(i) = camLocTable.MinApplicablePoint(i);
        maxAppPt(i) = camLocTable.MaxApplicablePoint(i);

        % Perform absolute orientation [Horn 1987] to find rotation and 
        % translation between total station and camera coordinate frame
        eigMat = zeros(4);
        for j = 1:size(tP_C,2)
            eigMat = eigMat + quatbarmat([tP_T(:,j); 0])' * ...
                        quatmat([tP_C(:,j); 0]);
        end
        [eigVec,eigVal] = eig(eigMat);
        [~,indMax]      = max(diag(eigVal));
        qTtoC           = quatconj(eigVec(:,indMax));
        RTtoC(:,:,i)    = quat2rot(qTtoC);
        tTtoC_C         = RTtoC(:,:,i) * tPm_T - tPm_C;
        tCtoT_C(:,i)    = -tTtoC_C;
        
        % Measure the error in each points' reprojection
        err     = tPo_C - (RTtoC(:,:,i)*tPo_T - tTtoC_C);
        errPt   = sqrt(err(1,:).^2 + err(2,:).^2 + err(3,:).^2);
        disp('Rotation RTtoC:');
        RTtoC(:,:,i)
        disp(['Norm of reprojection error for each point with '     ...
                'applicable total stn points ' num2str(minAppPt(i)) ...
                ' to ' num2str(maxAppPt(i)) ' [m]: ' num2str(errPt)]);
    end

    % Loop through images
    imDims = [];
    tR = []; tG = []; tB = [];      % filter transmissivities
    tGtoR = []; tBtoR = [];         % filter crosstalk to red img channel
    tRtoG = []; tBtoG = [];         % filter crosstalk to green img channel
    tRtoB = []; tGtoB = [];         % filter crosstalk to blue img channel
    pRc = []; pGc = []; pBc = [];   % coordinates of PSF centroids [px]
    DphR = []; DphG = []; DphB = [];% pinhole diameters [px]
    dispBG = []; dispBR = []; dispRG = [];  % disparities at each depth
    imagesUsed = {};                % images used in final calibration
    tPts_C = [];                    % total station points used
    depths = [];                    % depths used in final calibration
    imPSFsUsed = {};                % PSFs used in final calibration
    for pt = 1:(length(calibTable.Range))

        % Read and show image
        imageName = [calibTable.ImagePSF{pt} '.' ext];
        im = im2double(imread([folder imageName]));
        imDims = size(im);
        figure(1); clf; imshow(im); hold on;

        % Get the measured depth from the total station points
        pointId = calibTable.PointId(pt);
        indApp  = find(pointId >= minAppPt & pointId <= maxAppPt,1);
        tPt_T   = tsPts(tsPts(:,1) == pointId,2:4)';
        tPt_C   = tCtoT_C(:,indApp) + RTtoC(:,:,indApp) * tPt_T;
        depth   = abs(tPt_C(3));

        % Ask user to either skip the image or put a box around it
        depthApprox = calibTable.Range(pt);
        resp = questdlg(['Image ' imageName ' taken at a depth of ~'    ...
            num2str(depthApprox) 'm. Total station point id '           ...
            num2str(pointId) ' measured depth of ' num2str(depth) 'm '  ...
            'from CCD centroid.','Drag a box around PSF at POINT'       ...
            'LOCATION ' num2str(calibTable.PointLoc(pt))                ...
            ', or skip this point.'],'Process Point','Proceed','Skip',  ...
            'Proceed');
        if strcmp(resp,'Skip')
            continue;
        end

        % Get user-input box around PSF, and extract image section
        hPSF   = imrect(gca);     
        posPSF = hPSF.getPosition();
        imPSF  = im(round(posPSF(2)):round(posPSF(2)+posPSF(4)), ...
                    round(posPSF(1)):round(posPSF(1)+posPSF(3)),:);
        redPSF = imPSF(:,:,1); greenPSF = imPSF(:,:,2); bluePSF = imPSF(:,:,3);

        % Whether to algorithmically skip this point
        skipPoint = 0;

        % Loop through circles, processing them
        theta = 0:0.01:2*pi;
        figure(1); clf; imshow(imPSF); hold on;
        for j = 1:3

            % Extract the strongest circle from the normalized PSF image
            imPSFg = imPSF(:,:,j)/max(max(imPSF(:,:,j)));
            [centers,radii,metric] = imfindcircles(imPSFg, ...
                [3 round(max(size(imPSF))/4)],'ObjectPolarity','bright', ...
                'Sensitivity',0.95);
            radii = radii - 1.5;

            if isempty(centers)
                skipPoint = 1;
                disp(['No circle found in color channel ' num2str(j) ...
                        ', for image ' imageName ', point ' ...
                        num2str(calibTable.PointLoc(pt)) ', skipping point!']);
                break;
            else

                % Create a circular mask around the strongest detected circles
                colCirc  = round(centers(1,1) + cos(theta)*(radii(1)));
                rowCirc  = round(centers(1,2) + sin(theta)*(radii(1)));
                circMask = roipoly(imPSF,colCirc,rowCirc);

                % Get the mean red, green, and blue values inside the circle
                redVal   = mean(redPSF(circMask));
                greenVal = mean(greenPSF(circMask));
                blueVal  = mean(bluePSF(circMask));

                % Get mean intensity in area (i.e. iGinAr means green intensity 
                % in red circle area), as well as center positions
                if max([redVal greenVal blueVal]) == redVal
                    iRinAr = redVal; iGinAr = greenVal; iBinAr = blueVal;
                    pRcHyp = centers(1,:)'; DphRhyp = radii(1,:)'; 
                    col = 'r';
                elseif max([redVal greenVal blueVal]) == greenVal
                    iGinAg = greenVal; iRinAg = redVal; iBinAg = blueVal;
                    pGcHyp = centers(1,:)'; DphGhyp = radii(1,:)'; 
                    col = 'g';
                elseif max([redVal greenVal blueVal]) == blueVal
                    iBinAb = blueVal; iRinAb = redVal; iGinAb = greenVal;
                    pBcHyp = centers(1,:)'; DphBhyp = radii(1,:)'; 
                    col = 'b';
                end

                % Show the circle
                viscircles(centers(1,:),radii(1),'Color',col);
                plot(centers(1,1),centers(1,2),'x','Color','w');
            end

        end

        set(gcf, 'position', [100,100,1024,1024]);

        % If skipping the point, display a dialog and move onto the next pt
        if skipPoint
            resp = questdlg(['No circle found in color channel '       ...
                num2str(j) ', for image ' imageName ', point '         ...
                num2str(calibTable.PointLoc(pt)) ', skipping point!'], ...
                'Skipping Point','Ok','Ok');
            continue;
        end

        % Get the transmissivities, normalizing by the amount of light that
        % landed in the intended location, such that tR + tG + tB = 1
        tRhyp = iRinAr/(iRinAr+iGinAg+iBinAb);
        tGhyp = iGinAg/(iRinAr+iGinAg+iBinAb);
        tBhyp = iBinAb/(iRinAr+iGinAg+iBinAb);

        % Get the crosstalk between channels as ratio between light that
        % landed in an unintented location to that which landed as intended
        % (This method is deprecated, use calibrateTriangleColorSpace.m)
        tGtoRhyp = iRinAg/iGinAg;   tBtoRhyp = iRinAb/iBinAb;
        tRtoGhyp = iGinAr/iRinAr;   tBtoGhyp = iGinAb/iBinAb;
        tRtoBhyp = iBinAr/iRinAr;   tGtoBhyp = iBinAg/iGinAg;

        % Add results if they appear accurate
        drawnow;
        goodToGo = questdlg({'Commit the following calibration data?', ...
                    ['tR: ' num2str(tRhyp) ', tG: ' num2str(tGhyp)     ...
                    ', tB: ' num2str(tBhyp) ],['tGtoR: '               ...
                    num2str(tGtoRhyp) ', tBtoR: ' num2str(tBtoRhyp)],  ...
                    ['tRtoG: ' num2str(tRtoGhyp) ', tBtoG: '           ...
                    num2str(tBtoGhyp)],['tRtoB: ' num2str(tRtoBhyp)    ...
                    ', tGtoB: ' num2str(tGtoBhyp) ],['pRc: '           ...
                    num2str(pRcHyp')],['pGc: ' num2str(pGcHyp')],      ...
                    ['pBc: ' num2str(pBcHyp')]},                       ...
                    'Commit Data','Yes','No','Yes');

        if strcmp(goodToGo,'Yes')
            tPts_C(:,end+1) = tPt_C;
            depths(end+1) = depth;
            imagesUsed{end+1} = imageName;
            imPSFsUsed{end+1} = imPSF;
            tR(end+1) = tRhyp; tG(end+1) = tGhyp; tB(end+1) = tBhyp;
            tGtoR(end+1) = tGtoRhyp;        tBtoR(end+1) = tBtoRhyp;
            tRtoG(end+1) = tRtoGhyp;        tBtoG(end+1) = tBtoGhyp;
            tRtoB(end+1) = tRtoBhyp;        tGtoB(end+1) = tGtoBhyp;
            pRc(:,end+1) = round(posPSF(1:2))' + pRcHyp;          
            pGc(:,end+1) = round(posPSF(1:2))' + pGcHyp;
            pBc(:,end+1) = round(posPSF(1:2))' + pBcHyp;
            DphR(:,end+1) = DphRhyp; 
            DphG(:,end+1) = DphGhyp; 
            DphB(:,end+1) = DphBhyp;
            dispBG(end+1) = sign(pGcHyp(2)-pBcHyp(2))*norm(pGcHyp-pBcHyp);
            dispBR(end+1) = sign(pRcHyp(1)-pBcHyp(1))*norm(pRcHyp-pBcHyp);
            dispRG(end+1) = sign(pGcHyp(2)-pRcHyp(2))*norm(pGcHyp-pRcHyp);           
        end     

    end

    % Attempt to fit calibration parameters to depth vs. disparity for each
    % channel separately
    x0 = [1 -0.1 f 0.75]'; % [ a c zD b ]'
    Acon = [0 0 0 -1 ;
            0 0 0  0 ;
            0 0 0  0 ];    % constrain baselines to be >= 0
    bcon = [0 0 0 ];
    opt = optimoptions('fmincon','StepTolerance',1e-10,         ...
                    'MaxFunctionEvaluations',1e5,'MaxIterations',1e4);
    xOpt = zeros(length(x0),3); cost = zeros(3,1);
    for j = 1:3
        [xOpt(:,j),cost(j)] = fmincon(@(x)                      ...
            calibCost(x,tPts_C,pRc,pGc,pBc,sPx,j),x0,Acon,bcon, ...
            [],[],[],[],[],opt); 
    end

    % Attempt to fit calibration parameters to green image centroid and
    % effective focal length, and gap between front and rear principal plane
    y0 = [imDims(2)/2 imDims(1)/2 f 0]'; % [ cGx cGy fEff gEff ]'
    Acon = [-1  0  0  0 ;
             0 -1  0  0 ;
             0  0 -1  0 ;
             0  0  0 -1];    % constrain all to be >= 0
    bcon = [0 0 0 0];
    opt = optimoptions('fmincon','StepTolerance',1e-10,         ...
                    'MaxFunctionEvaluations',1e5,'MaxIterations',1e5);
    [yOpt,costY] = fmincon(@(y)                      ...
            calibCost2(y,tPts_C,pGc,sPx),y0,[],[], ...
            [],[],[],[],[],opt);
    cGx = yOpt(1); cGy = yOpt(2); fEff = yOpt(3); gEff = yOpt(4);

    % Note that calibration was completed
    calibCompleted = 1;

% Loading data ----------------------------------------------------------
elseif loadData == 1
    
    dirSplit = strsplit(folder,'/');
    load([folder camName '-' dirSplit{end-1} '.mat']);
    
end

% Model depth vs. disparity with optimal parameters
if size(xOpt,2) == 1    % when loading calibration 06.01.2018
    aBG = xOpt(1); aBR = xOpt(1); aRG = xOpt(1);
    cBG = xOpt(2); cBR = xOpt(2); cRG = xOpt(2);
    zDBG = xOpt(3); zDBR = xOpt(3); zDRG = xOpt(3);
    bBG = xOpt(4); bBR = xOpt(5); bRG = xOpt(6);
else                    % all newer calibrations
    aBG = xOpt(1,1); aBR = xOpt(1,2); aRG = xOpt(1,3);
    cBG = xOpt(2,1); cBR = xOpt(2,2); cRG = xOpt(2,3);
    zDBG = xOpt(3,1); zDBR = xOpt(3,2); zDRG = xOpt(3,3);
    bBG = xOpt(4,1); bBR = xOpt(4,2); bRG = xOpt(4,3);
end
dispBGmod  = dispFromDepthPx(depths,aBG,bBG,cBG,zDBG,sPx);
dispBRmod  = dispFromDepthPx(depths,aBR,bBR,cBR,zDBR,sPx);
dispRGmod  = dispFromDepthPx(depths,aRG,bRG,cRG,zDRG,sPx);
depthBGmod = depthFromDispPx(dispBG,aBG,bBG,cBG,zDBG,sPx);
depthBRmod = depthFromDispPx(dispBR,aBR,bBR,cBR,zDBR,sPx);
depthRGmod = depthFromDispPx(dispRG,aRG,bRG,cRG,zDRG,sPx);

% Model reprojection with optimized parameters
pGcMod = zeros(size(pGc));
for j = 1:size(tPts_C,2)
    pGcMod(:,j) = pinholeProject(tPts_C(:,j),fEff,cGx,cGy,sPx,gEff);
end

% Printing data ---------------------------------------------------------

% Print optimization results
disp('Model --------------------------------------------------');
disp('    disparityPx = b / sPx * (1 - a + zD/(depth + c)');
disp('    depth     = zD / (disparityPx*sPx/b - 1 + a) - c');
disp(['Parameters aBG, aBR, aRG : ' num2str(aBG) ', ' num2str(aBR) ', ' ...
    num2str(aRG) ]);
disp(['Parameters cBG, cBR, cRG : ' num2str(cBG) ', ' num2str(cBR) ', ' ...
    num2str(cRG) ]);
disp(['Parameters zDBG, zDBR, zDRG : ' num2str(zDBG) ', '               ...
    num2str(zDBR) ', ' num2str(zDRG) ]);
disp(['Parameters bBG, bBR, bRG : ' num2str(bBG) ', ' num2str(bBR) ', ' ...
    num2str(bRG) ]);
disp(['Calibration optim costs  : ' num2str(cost')]);

% Get mean values of each parameter
nD = length(depths);
tRm = mean(tR); tGm = mean(tG); tBm = mean(tB);
tGtoRm = mean(tGtoR); tBtoRm = mean(tBtoR);
tRtoGm = mean(tRtoG); tBtoGm = mean(tBtoG);
tRtoBm = mean(tRtoB); tGtoBm = mean(tGtoB);
pRcm = zeros(2,1); pGcm = zeros(2,1); pBcm = zeros(2,1);
indExcl = find(strcmp(imNameExcludeCentroids,imagesUsed),1);
nIncl   = indExcl-1;
for i = 1:indExcl-1
    % Get distances from the center of red, and normalize by line R-G
    pGcm = pGcm + (pGc(:,i)-pRc(:,i))/norm(pGc(:,i)-pRc(:,i))*1/nIncl;
    pBcm = pBcm + (pBc(:,i)-pRc(:,i))/norm(pGc(:,i)-pRc(:,i))*1/nIncl;    
end

disp('Excel --------------------------------------------------');
disp(['tR tG tB tGtoR tBtoR tRtoG tBtoG tRtoB tGtoB pRcX pRcY '       ...
        'pGcX pGcY pBcX pBcY aBG aBR aRG cBG cBR cRG zDBG zDBR zDRG ' ...
        'bBG bBR bRG fEff cGx cGy gEff = ']);
fprintf(['%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'...
    ' %f %f %f %f %f %f %f'],tRm,tGm,tBm,tGtoRm,tBtoRm,tRtoGm,tBtoGm, ...
    tRtoBm,tGtoBm,pRcm(1),pRcm(2),pGcm(1),pGcm(2),pBcm(1),pBcm(2),    ...
    aBG,aBR,aRG,cBG,cBR,cRG,zDBG,zDBR,zDRG,bBG,bBR,bRG,fEff,cGx,cGy,gEff);

disp('   ');

disp('C++ ----------------------------------------------------');
fprintf(['double tR = %f, tG = %f, tB = %f,\n       tGtoR = %f, '     ...
    'tBtoR = %f, tRtoG = %f,\n       tBtoG = %f, tRtoB = %f,'         ...
    ' tGtoB = %f;\ncv::Point2f pRc(%f,%f), pGc(%f,%f),\n           '  ...
    'pBc(%f,%f);\ndouble aBG = %f, aBR = %f, aRG = %f,'               ...
    '\n       cBG = %f, cBR = %f, cRG = %f,'                          ...
    '\n       zDBG = %f, zDBR = %f, zDRG = %f, sPx = %e,'             ...
    '\n       bBG = %f, bBR = %f, bRG = %f, fEff = %f,'               ...
    '\n       cGx = %f, cGy = %f, gEff = %f;\n'],tRm,tGm,tBm,tGtoRm,  ...
    tBtoRm,tRtoGm,tBtoGm,tRtoBm,tGtoBm,pRcm(1),pRcm(2),pGcm(1),       ...
    pGcm(2),pBcm(1),pBcm(2),aBG,aBR,aRG,cBG,cBR,cRG,zDBG,zDBR,zDRG,   ...
    sPx,bBG,bBR,bRG,fEff,cGx,cGy,gEff);

disp('   ');
disp('MATLAB -------------------------------------------------');
fprintf(['tR = %f; tG = %f; tB  = %f;\ntGtoR = %f; tBtoR = %f;'         ...
    '\ntRtoG = %f; tBtoG = %f;\ntRtoB = %f; tGtoB = %f;\n'              ...
    'pRc = [%f,%f]''; pGc = [%f,%f]''; pBc = [%f,%f]'';\n'              ...
    'aBG = %f; aBR = %f; aRG = %f;\n'                                   ...
    'cBG = %f; cBR = %f; cRG = %f;\n'                                   ...
    'zDBG = %f; zDBR = %f; zDRG = %f;\n'                                ...
    'bBG = %f; bBR = %f; bRG = %f;\n'                                   ...
    'fEff = %f; cGx = %f; cGy = %f; gEff = %f;\n'],                     ...                                ...
    tRm,tGm,tBm,tGtoRm,tBtoRm,tRtoGm,tBtoGm,tRtoBm,tGtoBm,pRcm(1),      ...
    pRcm(2),pGcm(1),pGcm(2),pBcm(1),pBcm(2),aBG,aBR,aRG,cBG,cBR,cRG,    ...
    zDBG,zDBR,zDRG,bBG,bBR,bRG,fEff,cGx,cGy,gEff);

% Display crosstalk correction; limit the total number of images to show 
% figure(1); clf;
% nCross = min(nD,7);
% for i = 1:nCross
%     R = imPSFsUsed{i}(:,:,1); G = imPSFsUsed{i}(:,:,2); 
%     B = imPSFsUsed{i}(:,:,3);
%     subplot(nCross,7,7*(i-1)+1); imshow(imPSFsUsed{i});
%     subplot(nCross,7,7*(i-1)+2); imshow(R);
%     subplot(nCross,7,7*(i-1)+3); imshow(G);
%     subplot(nCross,7,7*(i-1)+4); imshow(B);
%     subplot(nCross,7,7*(i-1)+5); imshow(R - tGtoRm*G - tBtoRm*B);
%     subplot(nCross,7,7*(i-1)+6); imshow(G - tRtoGm*R - tBtoGm*B);
%     subplot(nCross,7,7*(i-1)+7); imshow(B - tRtoBm*R - tGtoBm*G);    
% end

% Display depth vs disparity
figure(1); clf;
depthDispPlot = [depths' dispBG' dispBR' dispRG' ...
                    dispBGmod' dispBRmod' dispRGmod'];
depthDispPlot = sortrows(depthDispPlot);
plot(depthDispPlot(:,1),depthDispPlot(:,2),'x','Color',0.5*[0 1 1]); 
    hold on; grid on;
plot(depthDispPlot(:,1),depthDispPlot(:,3),'d','Color',0.5*[1 0 1]);  
plot(depthDispPlot(:,1),depthDispPlot(:,4),'^','Color',0.5*[1 1 0]);
plot(depthDispPlot(:,1),depthDispPlot(:,5),'-','Color',0.5*[0 1 1]);
plot(depthDispPlot(:,1),depthDispPlot(:,6),'-.','Color',0.5*[1 0 1]);  
plot(depthDispPlot(:,1),depthDispPlot(:,7),'--','Color',0.5*[1 1 0]);
title('Depth vs Disparity'); xlabel('Depth [m]'); 
ylabel('Disparity or Diameter [px]');
legend('Meas Disparity BG','Meas Disparity BR','Meas Disparity RG', ...
    'Modeled Disparity BG','Modeled Disparity BR',                  ...
    'Modeled Disparity RG');

% Display total station points from the perspective of the camera
figure(2); clf; ax = 0.25;
plot3(0,0,0,'or'); hold on; grid on; axis equal;    % cam origin
plot3([0 0],[0 ax],[0 0],'-r');                     % cam x
plot3([0 0],[0 0],-[0 ax],'-g');                     % cam y
plot3(-[0 ax],[0 0],[0 0],'-b');                    % cam z
plot3(-tPts_C(3,:),tPts_C(1,:),-tPts_C(2,:),'xk');  % pinhole position
for i = 1:nCamPoses
    plot3(-tCtoT_C(3,i),tCtoT_C(1,i),-tCtoT_C(2,i),'xb');
    plot3([-tCtoT_C(3,i) -tCtoT_C(3,i)-RTtoC(3,1,i)*ax], ...
          [ tCtoT_C(1,i)  tCtoT_C(1,i)+RTtoC(1,1,i)*ax], ...
          [-tCtoT_C(2,i) -tCtoT_C(2,i)-RTtoC(2,1,i)*ax],'--r');  % total station x
    plot3([-tCtoT_C(3,i) -tCtoT_C(3,i)-RTtoC(3,2,i)*ax], ...
          [ tCtoT_C(1,i)  tCtoT_C(1,i)+RTtoC(1,2,i)*ax], ...
          [-tCtoT_C(2,i) -tCtoT_C(2,i)-RTtoC(2,2,i)*ax],'--g');  % total station y
    plot3([-tCtoT_C(3,i) -tCtoT_C(3,i)-RTtoC(3,3,i)*ax], ...
          [ tCtoT_C(1,i)  tCtoT_C(1,i)+RTtoC(1,3,i)*ax], ...
          [-tCtoT_C(2,i) -tCtoT_C(2,i)-RTtoC(2,3,i)*ax],'--b');  % total station z
end
title('Test Geometry'); legend('CCD centroid','x_{cam}','y_{cam}', ...
    'z_{cam}','Pinhole positions','Total Station Origin','x_{ts}', ...
    'y_{ts}','z_{ts}');
xlabel('-^Cz [m]'); ylabel('^Cx [m]'); zlabel('-^Cy [m]');

% Save the data if it is new (not loaded)
if loadData == 0 && calibCompleted == 1 && saveData
    dirSplit = strsplit(folder,'/');
    save([folder camName '-' dirSplit{end-1} '.mat'],'tPts_C','depths', ...
            'imagesUsed','imPSFsUsed','tR','tG','tB','tGtoR','tBtoR',   ...
            'tRtoG','tBtoG','tRtoB','tGtoB','pRc','pGc','pBc',          ...
            'dispBG','dispBR','dispRG','imDims','zR','sPx','xOpt',      ...
            'cost','RTtoC','tCtoT_C','DphR','DphB','DphG','nCamPoses',  ...
            'fEff','cGx','cGy','gEff');
    disp(['Calibration written to: ' folder camName '-' ...
            dirSplit{end-1} '.mat']);
end

% Display modeling error in disparity in depth vs disparity
figure(3); clf;
errMod = depthDispPlot(:,5:7)-depthDispPlot(:,2:4);
plot(depthDispPlot(:,1),errMod(:,1),'x','Color',0.5*[0 1 1]);
    hold on; grid on;
plot(depthDispPlot(:,1),errMod(:,2),'d','Color',0.5*[1 0 1]);  
plot(depthDispPlot(:,1),errMod(:,3),'^','Color',0.5*[1 1 0]);
title('Error in Model'); xlabel('Depth [m]'); 
ylabel('Modeled Disparity - Measured Disparity [px]');
legend('Disparity BG','Disparity BR','Disparity RG');

% Display modeling error in depth in depth vs disparity
figure(4); clf;
errModD = [depths-depthBGmod; depths-depthBRmod; depths-depthRGmod];
plot(depths,errModD(1,:),'x','Color',0.5*[0 1 1]);
    hold on; grid on;
plot(depths,errModD(2,:),'d','Color',0.5*[1 0 1]);  
plot(depths,errModD(3,:),'^','Color',0.5*[1 1 0]);
ylim([-0.7 0.2]);
title('Error in Depth Model'); xlabel('Depth [m]'); 
ylabel('Error (Depth - Modeled Depth) [m]');
legend('BG','BR','RG');

% Display modeling error in projection
figure(5); clf;
plot(pGc(1,:),pGc(2,:),'ob'); hold on; grid on;
plot(pGcMod(1,:),pGcMod(2,:),'xr');
plot(imDims(2)/2,imDims(1)/2,'+k','MarkerSize',10);
plot(cGx,cGy,'+g','MarkerSize',10);
xlim([0 imDims(2)]); ylim([0 imDims(1)]);
axis equal; ax = gca; ax.YDir = 'reverse';
xlabel('Image x [px]'); ylabel('Image y [px]'); 
legend('Measured pGc','Projected pGc','Image Centroid', ...
    'Green Principal Point (c_{Gx}, c_{Gy})'); 
title('Total Station Measurement Projection');
