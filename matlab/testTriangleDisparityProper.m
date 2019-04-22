% This script tests the functionality of the TriangleDisparity class using
% diffraction-accurate PSFs to simulate an image at a fixed depth
% Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Read the original image file
imFileName  = '../data/triangle-disparity-test-images/trid1.png';
% imFileName  = '../data/triangle-disparity-test-images/trid2.png';
% imFileName  = '../data/triangle-disparity-test-images/trid3.png';
% imFileName  = '../data/triangle-disparity-test-images/trid4.png';
% imFileName  = '../data/triangle-disparity-test-images/trid5.png';
imOrig      = im2double(imread(imFileName));

disp('------------------------------------------------------------------');
disp(imFileName);
disp('------------------------------------------------------------------');

% Transmissivities of the filters
tR = 0.7; tG = 0.7; tB = 0.7;

% Crosstalk between filters and image
TCtoF = eye(3);

% Select whether to add shot noise to images
addShotNoise = 1;

% Camera properties (aperture diameter, focal length, pixel size)
cam     = struct();
cam.f   = 0.05;
cam.D   = cam.f/1.4;
cam.sPx = 4.5e-6;

% Specify the scene depth [m]
zO      = -3;

% Specify the reference focal depth [m] and baseline [m] and get detector 
% position and nominal disparity [px]
zR      = -3.5;
b       = cam.D * 0.75 * sqrt(3)/2;
zD      = zR*cam.f/(zR + cam.f);
dispNom = b/cam.sPx * (1 - zD/cam.f - zD./zO);

% Create masks
Nmask  = 1024;      % mask has dimensions Nmask x Nmask [px]
maskR  = getMask('tridispred',Nmask);
maskG  = getMask('tridispgreen',Nmask);
maskB  = getMask('tridispblue',Nmask);

% Centroid coordinates for red, green, and blue circles in PSF (x,y) [rAp]
pRc = [0.6495 -0.375]'; pGc = [0 0.75]'; pBc = [-0.6495 -0.375]';

% Create TriangleDisparity object and add focused image and masks
triangleDisparity = TriangleDisparity(tR,tG,tB,pRc,pGc,pBc,cam.D, ...
    cam.f,cam.sPx,TCtoF);
triangleDisparity.setFocusedImage(imOrig);
triangleDisparity.setMasks(maskR,maskB,maskG);

% Blur image and plot results
triangleDisparity.blurImageDepth(zR,zO,addShotNoise);
figure(2); clf; triangleDisparity.plotBlurredImages();

% Plot masks
figure(1); clf; triangleDisparity.plotMasks();

% Prepare rotated images and plot results
triangleDisparity.prepareRotatedImages();
figure(3); clf; triangleDisparity.plotRotatedImages();

% Calculate disparity map and plot results
winSize             = 15;
dispRange           = [0 16];
uniqueThresh        = 15;

[dispMap,dispMapBG,dispMapBR,dispMapRG] = ...
    triangleDisparity.calculateDisparity(winSize,dispRange,uniqueThresh);
figure(4); clf; 
triangleDisparity.plotDisparityResults(dispMap,dispMapBG,dispMapBR, ...
                                        dispMapRG,dispNom);

% Display stats on disparity finding
disp('Results ----------------------------------------------------------');
triangleDisparity.displayDisparityStats(dispMap,dispMapBG,dispMapBR, ...
                                        dispMapRG,dispNom);
