% This script loads the .mat file output by calibrateTriangeDisparity.m and
% determines the filter colorspace in RGB coordinates
% Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)

% Image extension
ext = 'JPG';

% Camera name
camName = 'dalsa5100';

% Image folder
folder = '../data/calib-glass50-led-08.24.2018/';
% folder = '../data/calib-glass85-led-08.24.2018/';

% Load the previous calibration file
dirSplit = strsplit(folder,'/');
load([folder camName '-' dirSplit{end-1} '.mat']);

% Loop through images of PSFs
iRinAr = []; iGinAr = []; iBinAr = [];  % intensities in red circle area
iRinAg = []; iGinAg = []; iBinAg = [];  % intensities in green circle area
iRinAb = []; iGinAb = []; iBinAb = [];  % intensities in blue circle area
for i = 1:length(imPSFsUsed)
       
    % Read image
    imPSF  = imPSFsUsed{i};
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
                iRinAr(end+1) = redVal; 
                iGinAr(end+1) = greenVal; 
                iBinAr(end+1) = blueVal;
                col = 'r';
            elseif max([redVal greenVal blueVal]) == greenVal
                iGinAg(end+1) = greenVal; 
                iRinAg(end+1) = redVal; 
                iBinAg(end+1) = blueVal;
                col = 'g';
            elseif max([redVal greenVal blueVal]) == blueVal
                iBinAb(end+1) = blueVal; 
                iRinAb(end+1) = redVal; 
                iGinAb(end+1) = greenVal;
                col = 'b';
            end

            % Show the circle
            viscircles(centers(1,:),radii(1),'Color',col);
            plot(centers(1,1),centers(1,2),'x','Color','w');
        end
        
    end
        
    drawnow;

end

% Get normalized vectors based on total intensity
iRinArN = iRinAr ./ sqrt( iRinAr.^2 + iGinAr.^2 + iBinAr.^2 );
iGinArN = iGinAr ./ sqrt( iRinAr.^2 + iGinAr.^2 + iBinAr.^2 );
iBinArN = iBinAr ./ sqrt( iRinAr.^2 + iGinAr.^2 + iBinAr.^2 );

iRinAgN = iRinAg ./ sqrt( iRinAg.^2 + iGinAg.^2 + iBinAg.^2 );
iGinAgN = iGinAg ./ sqrt( iRinAg.^2 + iGinAg.^2 + iBinAg.^2 );
iBinAgN = iBinAg ./ sqrt( iRinAg.^2 + iGinAg.^2 + iBinAg.^2 );

iRinAbN = iRinAb ./ sqrt( iRinAb.^2 + iGinAb.^2 + iBinAb.^2 );
iGinAbN = iGinAb ./ sqrt( iRinAb.^2 + iGinAb.^2 + iBinAb.^2 );
iBinAbN = iBinAb ./ sqrt( iRinAb.^2 + iGinAb.^2 + iBinAb.^2 );

vR = [mean(iRinArN) mean(iGinArN) mean(iBinArN)]';
vG = [mean(iRinAgN) mean(iGinAgN) mean(iBinAgN)]';
vB = [mean(iRinAbN) mean(iGinAbN) mean(iBinAbN)]';

% Get the matrix to transform from CCD colorspace to filter colorspace
A       = [vR vG vB];
TCtoF   = A\eye(3);     % inverse of A
disp('Transformation from CCD colorspace to filter colorspace TCtoF:');
TCtoF

% Display crosstalk correction; limit the total number of images to show
figure(1); clf;
nCross = min(size(imPSFsUsed,2),7);
for i = 1:nCross
    R = imPSFsUsed{i}(:,:,1); G = imPSFsUsed{i}(:,:,2); 
    B = imPSFsUsed{i}(:,:,3);
    rgbF = reshape((TCtoF*[R(:)'; G(:)'; B(:)'])', size(R,1),size(R,2),3);
    subplot(nCross,7,7*(i-1)+1); imshow(imPSFsUsed{i});
    subplot(nCross,7,7*(i-1)+2); imshow(R);
    subplot(nCross,7,7*(i-1)+3); imshow(G);
    subplot(nCross,7,7*(i-1)+4); imshow(B);
    subplot(nCross,7,7*(i-1)+5); imshow(rgbF(:,:,1));
    subplot(nCross,7,7*(i-1)+6); imshow(rgbF(:,:,2));
    subplot(nCross,7,7*(i-1)+7); imshow(rgbF(:,:,3));    
end

% Plot color vectors
figure(2); clf;

plot3([0 1],[0 0],[0 0],'-.r'); hold on; grid on;
plot3([0 0],[0 1],[0 0],'-.g');
plot3([0 0],[0 0],[0 1],'-.b');

for i = 1:length(iRinArN)
   plot3([0 iRinArN(i)],[0 iGinArN(i)],[0 iBinArN(i)],'xr'); 
end
plot3([0 vR(1)],[0 vR(2)],[0 vR(3)],'-r','LineWidth',2);

for i = 1:length(iRinAgN)
   plot3([0 iRinAgN(i)],[0 iGinAgN(i)],[0 iBinAgN(i)],'xg');
end
plot3([0 vG(1)],[0 vG(2)],[0 vG(3)],'-g','LineWidth',2);

for i = 1:length(iRinAbN)
   plot3([0 iRinAbN(i)],[0 iGinAbN(i)],[0 iBinAbN(i)],'xb');
end
plot3([0 vB(1)],[0 vB(2)],[0 vB(3)],'-b','LineWidth',2);

[xS,yS,zS] = sphere(50);
xS(xS < 0 | yS < 0 | zS < 0) = NaN;
yS(xS < 0 | yS < 0 | zS < 0) = NaN;
zS(xS < 0 | yS < 0 | zS < 0) = NaN;
surf(xS,yS,zS,'FaceAlpha',0,'EdgeAlpha',0.175);

axis equal;
title('Color space vectors and average');
xlabel('CCD Red'); ylabel('CCD Green'); zlabel('CCD Blue');
view(110,35);
