classdef TriangleDisparity < handle
%TRIANGLEDISPARITY Finds a depth map from color-coded aperture image
%   Generates and analyzes images formed using the triangular aperture mask
%   forming the PSF below:
%                O  tR*R            R: red filter
%                                   G: green filter
%                                   B: blue filter
%       tB*B  O      O  tG*G
% 
%   This mask features three circular pinholes which transmit red, green,
%   and blue light, arranged in a triangular pattern. 
%       tR: transmissivity of red pinhole
%       tG: transmissivity of green pinhole
%       tB: transmissivity of blue pinhole
%   
%   The disparity between the blue and green images, between the blue and
%   red images, and between the red and green images is calculated. To use
%   conventional algorithms and find horizontal disparity, the image
%   is rotated so that the disparity is in the horizontal direction for
%   each "left"-"right" comparison.
%   Author: Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)
    
    properties        
        
        % Transmissivity of each pinhole
        tR; tG; tB;
        
        % Crosstalk from filter to image
        TCtoF;
        
        % Pinhole centroid locations in units of aperture radius
        pRc; pGc; pBc;
        
        % Diffraction accurate masks and PSFs (proper)
        maskR; maskG; maskB; psfR; psfG; psfB;
        
        % Focused image
        imFoc; 
        
        % Blurred red, green, and blue images
        R; B; G;
        
        % Unit vectors between aperture centers in image coordinates
        pBtoG; pBtoR; pRtoG;
        
        % Rotation angles to facilitate horizontal disparity finding
        betaBG; betaBR; betaRG;
        
        % Rotated versions of the red, green, and blue images
        BrotBG; GrotBG; BrotBR; RrotBR; RrotRG; GrotRG;
        
        % Image region for recovery of images rotated back
        regRotBG; regRotBR; regRotRG;
  
        % Image and valid disparity area dims M = numRows, N = numCols
        M; N; Mv; Nv;
        
        % Number of photons expected per channel
        numPhotons;
        
        % Camera structure with aperture diameter D [m], focal length 
        % f [m], and pixel size sPx [m]
        cam;
        
        % Value for invalid disparity and invalid depth
        dispInv; depthInv;
         
    end
    
    methods
        
        %% Constructor
        function obj = TriangleDisparity(tR,tG,tB,pRc,pGc,pBc,D,f, ...
                                            sPx,TCtoF)
                              
            obj.tR = tR; obj.tG = tG; obj.tB = tB;
            obj.pRc = pRc;  obj.pGc = pGc;  obj.pBc = pBc;
            obj.pBtoG = (pGc-pBc)/norm(pGc-pBc); 
            obj.pBtoR = (pRc-pBc)/norm(pRc-pBc); 
            obj.pRtoG = (pGc-pRc)/norm(pGc-pRc);
            obj.TCtoF = TCtoF;
            
            % Get the angle to rotate the image counterclockwise by in deg
            obj.betaBG = atan2d(obj.pBtoG(2),obj.pBtoG(1));
            obj.betaBR = atan2d(obj.pBtoR(2),obj.pBtoR(1));
            obj.betaRG = atan2d(obj.pRtoG(2),obj.pRtoG(1));
                           
            obj.cam.D = D; obj.cam.f = f; obj.cam.sPx = sPx;
            obj.numPhotons = 7500;
            
        end
        
        %% Set focused image
        function setFocusedImage(obj,imFoc)
            obj.imFoc   = imFoc;
            obj.M       = size(imFoc,1);
            obj.N       = size(imFoc,2);
        end
        
        %% Set image (real image)
        function setImage(obj,im)
           
            % Image in color space
            Rc = im(:,:,1); Gc = im(:,:,2); Bc = im(:,:,3);
            
            % Transform to filter color space
            RGBc  = reshape( (obj.TCtoF*[Rc(:)'; Gc(:)'; Bc(:)'])',   ...
                                size(Rc,1),size(Rc,2),3 );
            obj.R = max(0,min(1, RGBc(:,:,1) ));
            obj.G = max(0,min(1, RGBc(:,:,2) ));
            obj.B = max(0,min(1, RGBc(:,:,3) ));
            obj.M = size(im,1);
            obj.N = size(im,2);
        end
        
        %% Set high resolution masks
        function setMasks(obj,maskR,maskG,maskB)           
            obj.maskR = maskR;
            obj.maskB = maskG;
            obj.maskG = maskB;            
        end   
        
        %% Create blurred images based on a focal plane reference depth [m]
        %  and object depth [m]. Add shot noise if addShotNoise==1 (do not 
        %  add shot noise by def)
        function blurImageDepth(obj,zR,zO,addShotNoise)
            
            % Wavelengths of each mask [m], and weight of each wavelength
            lR = (600:20:750)*1e-9; WlR = ones(1,length(lR));
            lG = (495:20:570)*1e-9; WlG = ones(1,length(lG));
            lB = (450:20:495)*1e-9; WlB = ones(1,length(lB));
            
            % Detector position
            zD = zR*obj.cam.f/(zR + obj.cam.f);
            
            % Get the PSFs as Npsf x Npsf px images with pixel size dx
            initProper;
            Npsf = max(size(obj.maskR));
            [obj.psfR, dxR] = getPsf(obj.cam.D,obj.cam.f,zO,zD,Npsf,lR, ...
                WlR,obj.maskR,0);
            [obj.psfG, dxG] = getPsf(obj.cam.D,obj.cam.f,zO,zD,Npsf,lG, ...
                WlG,obj.maskG,0);
            [obj.psfB, dxB] = getPsf(obj.cam.D,obj.cam.f,zO,zD,Npsf,lB, ...
                WlB,obj.maskB,0);
            
            % Resize the PSFs to match the camera pixel size
            obj.psfR = imresize(obj.psfR,dxR/obj.cam.sPx); 
            obj.psfG = imresize(obj.psfG,dxG/obj.cam.sPx); 
            obj.psfB = imresize(obj.psfB,dxB/obj.cam.sPx);
            
            % Normalize PSFs, but scale by maximum transmisivitty
            obj.psfR = obj.tR * obj.psfR/sum(obj.psfR(:)); 
            obj.psfG = obj.tG * obj.psfG/sum(obj.psfG(:));
            obj.psfB = obj.tB * obj.psfB/sum(obj.psfB(:));
            
            % Place the PSFs in a camera-sized image, assuming scaled PSF
            % is smaller than the image, and all PSFs are the same size            
            hPsf = size(obj.psfR,1);
            psfRs = zeros(obj.M,obj.N); psfGs = zeros(obj.M,obj.N);
            psfBs = zeros(obj.M,obj.N);
            ind1  = obj.M/2-floor(hPsf/2)+1:obj.M/2+ceil(hPsf/2);
            ind2  = obj.N/2-floor(hPsf/2)+1:obj.N/2+ceil(hPsf/2);
            psfRs(ind1,ind2) = obj.psfR; psfGs(ind1,ind2) = obj.psfG;
            psfBs(ind1,ind2) = obj.psfB;            
            
            % Get padded optical transfer functions (OTFs)
            otfR = fft2(fftshift(padarray(psfRs,[obj.M/2 obj.N/2])));
            otfG = fft2(fftshift(padarray(psfGs,[obj.M/2 obj.N/2])));
            otfB = fft2(fftshift(padarray(psfBs,[obj.M/2 obj.N/2])));
            
            % Get blurred image from convolution of image channels with OTF
            imPadR = ifft2( ...
                fft2(padarray(obj.imFoc(:,:,1),[obj.M/2 obj.N/2])).*otfR);
            imPadG = ifft2( ...
                fft2(padarray(obj.imFoc(:,:,2),[obj.M/2 obj.N/2])).*otfG);
            imPadB = ifft2( ...
                fft2(padarray(obj.imFoc(:,:,3),[obj.M/2 obj.N/2])).*otfB);
            
            % Remove image padding
            obj.R = imPadR(obj.M-obj.M/2+1:obj.M+obj.M/2, ...
                            obj.N-obj.N/2+1:obj.N+obj.N/2);
            obj.G = imPadG(obj.M-obj.M/2+1:obj.M+obj.M/2, ...
                            obj.N-obj.N/2+1:obj.N+obj.N/2);
            obj.B = imPadB(obj.M-obj.M/2+1:obj.M+obj.M/2, ...
                            obj.N-obj.N/2+1:obj.N+obj.N/2); 
                        
            % Add shot noise if requested
            if nargin > 3 && addShotNoise == 1
                disp(['Adding shot noise to image with ' ...
                    num2str(obj.numPhotons) ' photons/px/channel']);
                obj.R = poissrnd(max(0,min(1,obj.R))*obj.numPhotons) ...
                                            /obj.numPhotons;
                obj.B = poissrnd(max(0,min(1,obj.B))*obj.numPhotons) ...
                                            /obj.numPhotons;
                obj.G = poissrnd(max(0,min(1,obj.G))*obj.numPhotons) ...
                                            /obj.numPhotons;
            end
                        
        end
        
        %% Prepare the rotated images for future use
        function prepareRotatedImages(obj)
            
            tic;
                                  
            % Rotate the image using bilinear interpolation
            if abs(obj.betaBG) > 1e-2
                obj.BrotBG = imrotate(obj.B,obj.betaBG,'loose','bilinear');
                obj.GrotBG = imrotate(obj.G,obj.betaBG,'loose','bilinear');
                obj.regRotBG = imrotate(imrotate(ones(size(obj.B)), ...
                    obj.betaBG,'loose','bilinear'),-obj.betaBG);
            else
                obj.BrotBG   = obj.B;
                obj.GrotBG   = obj.G;
                obj.regRotBG = ones(size(obj.B));
            end
            if abs(obj.betaBR) > 1e-2
                obj.BrotBR = imrotate(obj.B,obj.betaBR,'loose','bilinear');
                obj.RrotBR = imrotate(obj.R,obj.betaBR,'loose','bilinear');
                obj.regRotBR = imrotate(imrotate(ones(size(obj.B)), ...
                    obj.betaBR),-obj.betaBR);
            else
                obj.BrotBR   = obj.B;
                obj.RrotBR   = obj.R;
                obj.regRotBR = ones(size(obj.B));
            end
            if abs(obj.betaRG) > 1e-2
                obj.RrotRG = imrotate(obj.R,obj.betaRG,'loose','bilinear');
                obj.GrotRG = imrotate(obj.G,obj.betaRG,'loose','bilinear');
                obj.regRotRG = imrotate(imrotate(ones(size(obj.R)), ...
                    obj.betaRG),-obj.betaRG);
            else
                obj.RrotRG   = obj.R;
                obj.GrotRG   = obj.G;
                obj.regRotRG = ones(size(obj.R));
            end
               
            disp(['TriangleDisparity::prepareRotatedImages() ' ...
                'completed in ' num2str(toc) ' s']);
            
        end 
        
        %% Calculate the disparity map with s-by-s px window (s is odd)
        %  disparity range dRange = [dMin dMax], and uniqueness threshold.
        function [dispMap,dispMapBG,dispMapBR,dispMapRG,validMap] = ...
            calculateDisparity(obj,s,dispRange,uniqueThresh)
           
            tic;
            
            % Set invalid disparity to one less than the min disparity
            obj.dispInv = dispRange(1) - 1;
            
            % Choose the disparity method ('SemiGlobal' or 'BlockMatching')
            % Note that for some reason the max disparity range for block
            % matching is 16
            if dispRange(2) - dispRange(1) <= 16
                dispMethod = 'BlockMatching';
            else
                dispMethod = 'SemiGlobal';
            end
            
            % Calculate disparity between rotated blue and green channels
            dispMapBG = disparity(obj.GrotBG,obj.BrotBG,'Method',    ...
                dispMethod,'DisparityRange',dispRange,'BlockSize',s, ...
                'UniquenessThreshold',uniqueThresh);
            dispMapBR = disparity(obj.RrotBR,obj.BrotBR,'Method',    ...
                dispMethod,'DisparityRange',dispRange,'BlockSize',s, ...
                'UniquenessThreshold',uniqueThresh);
            dispMapRG = disparity(obj.GrotRG,obj.RrotRG,'Method',    ...
                dispMethod,'DisparityRange',dispRange,'BlockSize',s, ...
                'UniquenessThreshold',uniqueThresh);
            dispMapBG = double(dispMapBG);
            dispMapBR = double(dispMapBR);
            dispMapRG = double(dispMapRG);
            
            % Rotate disparity maps back to upright using the nearest 
            % neighbour method, so as not to potentially average valid and
            % invalid disparities
            dispMapBG = imrotate(dispMapBG,-obj.betaBG);
            dispMapBR = imrotate(dispMapBR,-obj.betaBR);
            dispMapRG = imrotate(dispMapRG,-obj.betaRG);
            dispMapBG = dispMapBG(any(obj.regRotBG'),any(obj.regRotBG));
            dispMapBR = dispMapBR(any(obj.regRotBR'),any(obj.regRotBR));
            dispMapRG = dispMapRG(any(obj.regRotRG'),any(obj.regRotRG));
            
            % If the image rotation creates an extra pixel, truncate
            dispMapBG = dispMapBG(1:obj.M,1:obj.N);
            dispMapBR = dispMapBR(1:obj.M,1:obj.N);
            dispMapRG = dispMapRG(1:obj.M,1:obj.N);
            
            % Change invalid disparities to dispInv, not -realmax('single')
            dispMapBG(dispMapBG == -realmax('single')) = obj.dispInv;
            dispMapBR(dispMapBR == -realmax('single')) = obj.dispInv;
            dispMapRG(dispMapRG == -realmax('single')) = obj.dispInv;
            
            % Make the highest disparity and zero disparity invalid, since 
            % they are often the source of errors
            dispMapBG((dispMapBG>=dispRange(2)-1) | (dispMapBG==0)) = ...
                                                            obj.dispInv;
            dispMapBR((dispMapBR>=dispRange(2)-1) | (dispMapBR==0)) = ...
                                                            obj.dispInv;
            dispMapRG((dispMapRG>=dispRange(2)-1) | (dispMapRG==0)) = ...
                                                            obj.dispInv;
            
            % Get valid values
            isValidBG = dispMapBG ~= obj.dispInv;
            isValidBR = dispMapBR ~= obj.dispInv;
            isValidRG = dispMapRG ~= obj.dispInv;
            
            % Get shifted dispMapBR results, moving them from red pinhole 
            % coords (xR,yR) to green pinhole coords (xG,yG)
            dispMapBRs = ones(obj.M,obj.N)*obj.dispInv; % shifted dispMapBR
            nAvgBRs    = zeros(obj.M,obj.N);            % num of avg'd vals
            for yR = 1:size(dispMapBR,1)
                for xR = 1:size(dispMapBR,2)
                    
                    % If disparity is valid, shift it from pinhole R coords
                    % to pinhole G by the value of disparity
                    if isValidBR(yR,xR)
                        xG = round(xR + dispMapBR(yR,xR)*obj.pRtoG(1));
                        yG = round(yR + dispMapBR(yR,xR)*obj.pRtoG(2));
                        
                        % Check new location is in image, and if so, add
                        % the disparity to the shifted result
                        if xG >= 1 && xG <= obj.N && yG >= 1 && yG <= obj.M
                            nAvgBRs(yG,xG) = nAvgBRs(yG,xG) + 1;
                            if nAvgBRs(yG,xG) == 1
                                dispMapBRs(yG,xG) = dispMapBR(yR,xR);
                            else
                                dispMapBRs(yG,xG) = dispMapBRs(yG,xG)  ...
                                  *(nAvgBRs(yG,xG)-1)/nAvgBRs(yG,xG)   ...
                                  + dispMapBR(yR,xR)/nAvgBRs(yG,xG)   ;                                                   
                            end
                        end                        
                    end                    
                    
                end
            end
            
            % Reassign dispMapBR to shifted version; re-check validity
            dispMapBR = dispMapBRs;
            isValidBR = dispMapBR ~= obj.dispInv;
            
            % Create an M x N x 3 validity image, which is a three layer 
            % logical matrix, with a value of 0 when the disparity map is
            % invalid at that pixel, and 1 when it is valid. The three
            % layers refer to validity of dispMapBG, dispMapBR, and
            % dispMapRG, respectively.
            validMap = zeros(obj.M,obj.N,3);
                      
            % Determine which disparity maps contain valid results
            validMap(:,:,1) = isValidBG; % BG
            validMap(:,:,2) = isValidBR; % BR
            validMap(:,:,3) = isValidRG; % RG
            
            % Determine the number of valid disparities at each pixel
            numValid = sum(validMap,3);
            
            % Average the valid disparities
            dispMap  = (dispMapBG .* validMap(:,:,1) ...
                      + dispMapBR .* validMap(:,:,2) ...
                      + dispMapRG .* validMap(:,:,3) ) ./ numValid;
            dispMap(dispMap == Inf | dispMap == -Inf ...
                        | isnan(dispMap)) = obj.dispInv;

            figure(21); imshow(validMap); title('validMap');
                                    
            disp(['TriangleDisparity::calculateDisparity() ' ...
                'completed in ' num2str(toc) ' s']);
            
        end
        
        %% Calculate the depth given disparity, validity map and 
        %  calibration parameters
        function depthMap = calculateDepth(obj,dispMapBG,dispMapBR, ...
            dispMapRG,validMap,aBG,aBR,aRG,cBG,cBR,cRG,zDBG,zDBR,   ...
            zDRG,bBG,bBR,bRG)
           
            % Determine the number of valid disparities at each pixel
            numValid = sum(validMap,3);
            
            % Set a value for invalid depth
            obj.depthInv = -1;
            
            % Average the depth from valid disparities
            % depth = zD / (disparityPx*sPx/b - 1 + a) - c
            sPx = obj.cam.sPx;
            depthMap = ...
              ( (zDBG./(dispMapBG*sPx/bBG-1+aBG)-cBG).*validMap(:,:,1)  ...
              + (zDBR./(dispMapBR*sPx/bBR-1+aBR)-cBR).*validMap(:,:,2)  ...
              + (zDRG./(dispMapRG*sPx/bRG-1+aRG)-cRG).*validMap(:,:,3) )...
                ./ numValid;
            depthMap(depthMap == Inf | depthMap == -Inf ...
                        | isnan(depthMap)) = obj.depthInv;
                        
        end
        
        %% Calculate the point cloud x, y, and z are (MxN)
        function [ptCloud] = calculatePointCloud(obj,depthMap,fEff, ...
                            cGx,cGy,gEff,sPx)
           
            Md = size(depthMap,1); Nd = size(depthMap,2);
            u = repmat(1:Nd,Md,1); v = repmat((1:Md)',1,Nd); 
            
            ptCloud(:,:,1) = (depthMap-gEff) .* (u-cGx)*sPx/fEff;
            ptCloud(:,:,2) = (depthMap-gEff) .* (v-cGy)*sPx/fEff;
            ptCloud(:,:,3) = depthMap;                
            
        end
        
        %% Get the blurred image
        function imBlur = getBlurredImage(obj)
           imBlur(:,:,1) = obj.R;
           imBlur(:,:,2) = obj.G;
           imBlur(:,:,3) = obj.B;
        end
        
        %% Plot the masks at the full resolution provided
        function plotMasks(obj)
            im(:,:,1) = obj.maskR; im(:,:,2) = obj.maskG; 
            im(:,:,3) = obj.maskB;
            p(:,:,1) = obj.psfR/max([obj.psfR(:); obj.psfG(:); obj.psfB(:)]);
            p(:,:,2) = obj.psfG/max([obj.psfR(:); obj.psfG(:); obj.psfB(:)]);
            p(:,:,3) = obj.psfB/max([obj.psfR(:); obj.psfG(:); obj.psfB(:)]);
            subplot(2,4,1); imshow(im);        title('mask');
            subplot(2,4,2); imshow(obj.maskR); title('maskR');
            subplot(2,4,3); imshow(obj.maskG); title('maskG');
            subplot(2,4,4); imshow(obj.maskB); title('maskB');
            subplot(2,4,5); imshow(p);         title('PSF');
            subplot(2,4,6); imagesc(obj.psfR); 
                title('psfR'); axis equal;
            subplot(2,4,7); imagesc(obj.psfG);
                title('psfG'); axis equal;
            subplot(2,4,8); imagesc(obj.psfB);
                title('psfB'); axis equal;
        end
        
        %% Plot the blurred images
        function plotBlurredImages(obj)
            subplot(2,4,1); imshow(obj.imFoc); title('Focused');
            subplot(2,4,2); imshow(obj.R);     title('R Blurred');
            subplot(2,4,3); imshow(obj.B);     title('B Blurred');
            subplot(2,4,4); imshow(obj.G);     title('G Blurred');
            subplot(2,4,5); obj.plotBlurredImage();
        end
        
        %% Plot the blurred image
        function plotBlurredImage(obj)
           imBlur = obj.getBlurredImage();
           imshow(imBlur); title('Blurred Image');
        end
        
        %% Plot the scaled images
        function plotRotatedImages(obj)
            subplot(3,2,1); imshow(obj.BrotBG);     title('BrotBG');
            subplot(3,2,2); imshow(obj.GrotBG);     title('GrotBG');
            subplot(3,2,3); imshow(obj.BrotBR);     title('BrotBR');
            subplot(3,2,4); imshow(obj.RrotBR);     title('RrotBR');
            subplot(3,2,5); imshow(obj.RrotRG);     title('RrotRG');
            subplot(3,2,6); imshow(obj.GrotRG);     title('GrotRG');
        end
            
        %% Plot the disparity results 
        function plotDisparityResults(obj,dispMap,dispMapBG,dispMapBR, ...
                                        dispMapRG,dispNom)
            
            if nargin > 5
                dispErr   = abs(dispMap-dispNom);
                dispErrBG = abs(dispMapBG-dispNom);
                dispErrBR = abs(dispMapBR-dispNom);
                dispErrRG = abs(dispMapRG-dispNom);
                maxErr  = max([dispErr(:); dispErrBG(:);  ...
                                dispErrBR(:); dispErrRG(:)]);
            end
            
            minDisp = min([dispMap(:); dispMapBG(:);  ...
                            dispMapBR(:); dispMapRG(:)]);
            maxDisp = max([dispMap(:); dispMapBG(:);  ...
                            dispMapBR(:); dispMapRG(:)]);
            subplot(2,4,1); imagesc(dispMap); 
                title('Combined Disparity'); axis equal; colorbar;
                caxis([minDisp maxDisp]);
            subplot(2,4,2); imagesc(dispMapBG); title('Disparity Map BG');
                axis equal; colorbar; caxis([minDisp maxDisp]);
            subplot(2,4,3); imagesc(dispMapBR); title('Disparity Map BR');
                axis equal; colorbar; caxis([minDisp maxDisp]);
            subplot(2,4,4); imagesc(dispMapRG); title('Disparity Map RG');
                axis equal; colorbar; caxis([minDisp maxDisp]);
                
            if nargin > 5
                subplot(2,4,5); imagesc(dispErr); 
                    title('Combined Disparity Error');
                    axis equal; colorbar; caxis([0 maxErr]);
                subplot(2,4,6); imagesc(dispErrBG); 
                    title('Disparity BG Error');
                    axis equal; colorbar; caxis([0 maxErr]);
                subplot(2,4,7); imagesc(dispErrBR); 
                    title('Disparity BR Error');
                    axis equal; colorbar; caxis([0 maxErr]);
                subplot(2,4,8); imagesc(dispErrRG); 
                    title('Disparity RG Error');
                    axis equal; colorbar; caxis([0 maxErr]);
            end
        end
        
        %% Plot the depth results
        function plotDepthResults(obj,depthMap,minDepth,maxDepth)
            imagesc(depthMap); axis equal; colorbar; title('Depth Map');
                caxis([minDepth maxDepth]);
        end
        
        %% Plot a point cloud with depth results
        function plotPointCloud(obj,ptCloud,minDepth,maxDepth)
            z = ptCloud(:,:,3);
            z(z < minDepth | z > maxDepth) = obj.depthInv;
            ptCloud(:,:,3) = z;
            imDisp(:,:,1) = obj.G; imDisp(:,:,2) = obj.G; 
            imDisp(:,:,3) = obj.G; imDisp = histeq(imDisp);
            pcshow(ptCloud,imDisp); grid on; axis equal; axis off;
            title('Point Cloud');
            sW = obj.N*obj.cam.sPx; sH = obj.M*obj.cam.sPx; % sensor w & h            
            zlim([minDepth maxDepth]);
            xlim([-maxDepth*sW/2/obj.cam.f maxDepth*sW/2/obj.cam.f]);
            ylim([-maxDepth*sH/2/obj.cam.f maxDepth*sH/2/obj.cam.f]);  
            xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        end
        
        %% Display stats on disparity
        function [pcValidPx,avgDisp,avgDispErr] = displayDisparityStats(...
                        obj,dispMap,dispMapBG,dispMapBR,dispMapRG,dispNom)
                    
            dispErr         = abs(dispMap-dispNom);
            dispErrBG       = abs(dispMapBG-dispNom);
            dispErrBR       = abs(dispMapBR-dispNom);
            dispErrRG       = abs(dispMapRG-dispNom);
            
            validPx         = dispMap ~= obj.dispInv;
            numValidPx      = sum(validPx(:));
            pcValidPx       = numValidPx/(obj.M*obj.N) * 100;
            dispErrValidPx  = dispErr(validPx);
            dispMapValidPx  = dispMap(validPx);
            avgDisp         = mean(dispMapValidPx(:));
            avgDispErr      = sum(dispErrValidPx(:))/numValidPx;            
            disp(['Number of valid pixels                      : ' ...
                    num2str(numValidPx) '(' num2str(pcValidPx) '%)']);
            disp(['Average disparity for valid px              : ' ...
                    num2str(avgDisp)]);
            disp(['Disparity std dev for valid px              : ' ...
                    num2str(std(dispMapValidPx(:)))]);
            disp(['Average disparity error for valid pixels    : ' ...
                num2str(avgDispErr)]);
            disp(['Disparity error std dev for valid pixels    : ' ...
                num2str(std(dispErrValidPx(:)))]);
            
            validPx         = dispMapBG ~= obj.dispInv;
            numValidPx      = sum(validPx(:));
            pcValidPx       = numValidPx/(obj.M*obj.N) * 100;
            dispErrValidPx  = dispErrBG(validPx);
            dispMapValidPx  = dispMapBG(validPx);
            avgDisp         = mean(dispMapValidPx(:));
            avgDispErr      = sum(dispErrValidPx(:))/numValidPx;
            disp(['Number of valid pixels BG                   : ' ...
                    num2str(numValidPx) '(' num2str(pcValidPx) '%)']);
            disp(['Average disparity for valid px BG           : ' ...
                    num2str(avgDisp)]);
            disp(['Disparity std dev for valid px BG           : ' ...
                    num2str(std(dispMapValidPx(:)))]);
            disp(['Average disparity error for valid pixels BG : ' ...
                num2str(avgDispErr)]);
            disp(['Disparity error std dev for valid pixels BG : ' ...
                num2str(std(dispErrValidPx(:)))]);
                        
            validPx         = dispMapBR ~= obj.dispInv;
            numValidPx      = sum(validPx(:));
            pcValidPx       = numValidPx/(obj.M*obj.N) * 100;
            dispErrValidPx  = dispErrBR(validPx);
            dispMapValidPx  = dispMapBR(validPx);
            avgDisp         = mean(dispMapValidPx(:));
            avgDispErr      = sum(dispErrValidPx(:))/numValidPx;
            disp(['Number of valid pixels BR                   : ' ...
                    num2str(numValidPx) '(' num2str(pcValidPx) '%)']);
            disp(['Average disparity for valid px BR           : ' ...
                    num2str(avgDisp)]);
            disp(['Disparity std dev for valid px BR           : ' ...
                    num2str(std(dispMapValidPx(:)))]);
            disp(['Average disparity error for valid pixels BR : ' ...
                num2str(avgDispErr)]);
            disp(['Disparity error std dev for valid pixels BR : ' ...
                num2str(std(dispErrValidPx(:)))]);
            
            validPx         = dispMapRG ~= obj.dispInv;
            numValidPx      = sum(validPx(:));
            pcValidPx       = numValidPx/(obj.M*obj.N) * 100;
            dispErrValidPx  = dispErrRG(validPx);
            dispMapValidPx  = dispMapRG(validPx);
            avgDisp         = mean(dispMapValidPx(:));
            avgDispErr      = sum(dispErrValidPx(:))/numValidPx;
            disp(['Number of valid pixels RG                   : ' ...
                    num2str(numValidPx) '(' num2str(pcValidPx) '%)']);
            disp(['Average disparity for valid px RG           : ' ...
                    num2str(avgDisp)]);
            disp(['Disparity std dev for valid px RG           : ' ...
                    num2str(std(dispMapValidPx(:)))]);
            disp(['Average disparity error for valid pixels RG : ' ...
                num2str(avgDispErr)]);
            disp(['Disparity error std dev for valid pixels RG : ' ...
                num2str(std(dispErrValidPx(:)))]);
            
        end
    end
end

