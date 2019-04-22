% This script tests the functionality of the TriangleDisparity class using
% real images taken using a color-coded aperture camera
% Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Read the original image file
imFileName = '../data/blocks-glass50-09.26.2018/Dalsa-18.jpg';
% imFileName = '../data/blocks-glass85-09.26.2018/Dalsa-22.jpg';
im          = im2double(imread(imFileName));

disp('------------------------------------------------------------------');
disp(imFileName);
disp('------------------------------------------------------------------');

% Resize if desired
resizeFactor    = 1/4;
if resizeFactor ~= 1
    im = imresize(im,resizeFactor);
end

% Aperture calibration parameters 
if ~contains(imFileName,'glass85')
    % Aperture Glass50 - Calibration from 08.24.2018
    disp('Using glass50 calibration');
    tR = 0.368918; tG = 0.244706; tB  = 0.386376;
    tGtoR = 0.255327; tBtoR = 0.095088;
    tRtoG = 0.431083; tBtoG = 0.259603;
    tRtoB = 0.266450; tGtoB = 0.466142;
    pRc = [0.000000,0.000000]'; pGc = [-0.456487,0.889677]'; 
    pBc = [-1.026158,0.054581]';
    aBG = 1.001148; aBR = 1.000369; aRG = 1.001049;
    cBG = -0.063799; cBR = -0.064970; cRG = -0.060415;
    zDBG = 0.003363; zDBR = 0.001110; zDRG = 0.003052;
    bBG = 0.293332; bBR = 0.893171; bRG = 0.323857;
    fEff = 0.054958; cGx = 2456.156591; cGy = 2409.236163; 
    gEff = -0.020838;
    TCtoF = [ 1.259430503189212  -0.296170305552113  -0.042522437167015 ;
             -0.506915825715605   1.402281343425590  -0.305889481616210 ;
             -0.083531675722062  -0.523441191930818   1.183625763845906 ];
    minDepth = 1.5; maxDepth = 3; % min and max depth for Dalsa-18.jpg
else
    % Aperture Glass85 - Calibration from 08.24.2018
    disp('Using glass85 calibration');
    tR = 0.358558; tG = 0.289431; tB  = 0.352011;
    tGtoR = 0.261476; tBtoR = 0.098389;
    tRtoG = 0.435794; tBtoG = 0.282516;
    tRtoB = 0.234029; tGtoB = 0.459022;
    pRc = [0.000000,0.000000]'; pGc = [-0.760614,0.648781]'; 
    pBc = [-0.979609,-0.310073]';
    aBG = 1.000389; aBR = 1.000205; aRG = 1.000229;
    cBG = -0.048268; cBR = -0.059000; cRG = -0.059806;
    zDBG = 0.003177; zDBR = 0.001761; zDRG = 0.001910;
    bBG = 0.735652; bBR = 1.310770; bRG = 1.207753;
    fEff = 0.095356; cGx = 1863.453699; cGy = 3102.308784; 
    gEff = -0.210120;
    TCtoF = [ 1.257036294580806 -0.302448707745906 -0.035049432869650 ;
             -0.532219148021015  1.422040705546023 -0.339437934835568 ;
             -0.048405772799910 -0.529890758982390  1.200347639643513 ];
    minDepth = 3; maxDepth = 4.5; % min and max depth for Dalsa-22.jpg
end

% Modify calibration for rescaling
cGx = cGx * resizeFactor;
cGy = cGy * resizeFactor;

% Camera properties (aperture diameter, focal length, pixel size)
cam     = struct();
cam.f   = 0.05;
cam.D   = cam.f/1.4;
cam.sPx = 4.5e-6 / resizeFactor;

% Create TriangleDisparity object and add focused image and masks
triangleDisparity = TriangleDisparity(tR,tG,tB,pRc,pGc,pBc,cam.D, ...
    cam.f,cam.sPx,TCtoF);
triangleDisparity.setImage(im);

% Plot blurred images
figure(1); clf; triangleDisparity.plotBlurredImages();

% Prepare rotated images and plot results
triangleDisparity.prepareRotatedImages();
figure(2); clf; triangleDisparity.plotRotatedImages();

% Calculate disparity map and plot results
winSize             = 15;
dispRange           = [ -(ceil(resizeFactor*10)+1)*16 ...
                         (ceil(resizeFactor*10)+1)*16 ];
uniqueThresh        = 15;

[dispMap,dispMapBG,dispMapBR,dispMapRG,validMap] = ...
    triangleDisparity.calculateDisparity(winSize,dispRange,uniqueThresh);
figure(3); clf; 
triangleDisparity.plotDisparityResults(dispMap,dispMapBG,dispMapBR, ...
                                        dispMapRG);
 
% Calculate depth and plot results
depthMap = triangleDisparity.calculateDepth(dispMapBG,dispMapBR,       ...
            dispMapRG,validMap,aBG,aBR,aRG,cBG,cBR,cRG,zDBG,zDBR,zDRG, ...
            bBG,bBR,bRG);
figure(4); clf;
triangleDisparity.plotDepthResults(depthMap,minDepth,maxDepth);

% Get point cloud results and plot them
ptCloud = triangleDisparity.calculatePointCloud(depthMap,fEff, ...
                            cGx,cGy,gEff,cam.sPx);
figure(5); clf; 
triangleDisparity.plotPointCloud(ptCloud,minDepth,maxDepth);
