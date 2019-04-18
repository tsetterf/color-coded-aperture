/**
 * @file TriangleDisparity.cpp
 * @brief Finds a dense depth map from color-coded aperture image
 * @date May 18, 2018
 * @author Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)
 */

#include <fstream>
#include <float.h>
#include <math.h>
#include <algorithm>
#include <limits>
#include <ctime>
#include <cmath>
#include "TriangleDisparity.h"

/* ************************************************************************* */
TriangleDisparity::TriangleDisparity(double tR, double tG, double tB,
  cv::Mat TCtoF, cv::Point2f pRc, cv::Point2f pGc, cv::Point2f pBc, 
  int minDisparity, int numDisparities, int uniqueThresh, int windowSize)
  : tR_(tR), tG_(tG), tB_(tB), TCtoF_(TCtoF), pRc_(pRc), pGc_(pGc), pBc_(pBc),
  minDisparity_(minDisparity), maxDisparity_(minDisparity+numDisparities-1),
  numDisparities_(numDisparities), uniqueThresh_(uniqueThresh),
  windowSize_(windowSize), isCamCalibSet_(false), M_(0), N_(0), prepTime_(0.0),
  calcTime_(0.0) { 

  pBtoG_  = (pGc_-pBc_)/cv::norm(pGc_-pBc_);
  pBtoR_  = (pRc_-pBc_)/cv::norm(pRc_-pBc_);
  pRtoG_  = (pGc_-pRc_)/cv::norm(pGc_-pRc_);

  betaBG_ = std::atan2( pGc_.y-pBc_.y, pGc_.x-pBc_.x )*180.0/M_PI;
  betaBR_ = std::atan2( pRc_.y-pBc_.y, pRc_.x-pBc_.x )*180.0/M_PI;
  betaRG_ = std::atan2( pGc_.y-pRc_.y, pGc_.x-pRc_.x )*180.0/M_PI;
}


/* ************************************************************************* */
TriangleDisparity::~TriangleDisparity() { 
}

/* ************************************************************************* */
void TriangleDisparity::setImage(const cv::Mat& image) { 

  clock_t begin = clock(), end;

  // Convert channels to doubles (for crosstalk subtraction) if not already
  cv::Mat imageD;
  if (image.type() != CV_64FC3) { 
    image.convertTo(imageD,CV_64FC3,1.0);
  }

  // Convert image from CCD colorspace to filter colorspace
  cv::Mat imageF(imageD.rows,imageD.cols,CV_64FC3);
  for (int y = 0; y < image.rows; y++) { 

    const cv::Vec3d *imageDptr  = imageD.ptr<cv::Vec3d>(y);
    cv::Vec3d *imageFptr        = imageF.ptr<cv::Vec3d>(y);

    for (int x = 0; x < image.cols; x++) { 
      imageFptr[x] = cv::Vec3d( cv::Mat(TCtoF_ * cv::Mat(imageDptr[x])) );
    }

  }

  cv::Mat bgr[3];
  cv::split(imageF,bgr);

  // Convert the image to unsigned char, with appropriate intensity scaling
  bgr[2].convertTo(R_,CV_8UC1,1.0); 
  bgr[1].convertTo(G_,CV_8UC1,1.0);
  bgr[0].convertTo(B_,CV_8UC1,1.0);
  
  // If the image size has changed, create the bounding rotation rectangles
  if (image.rows != M_ || image.cols != N_) { 

    cv::Point2f center(0.5*((double) image.cols-1),0.5*((double) image.rows-1));

    regRotBG_ = cv::RotatedRect(center,image.size(),betaBG_);
    regRotBR_ = cv::RotatedRect(center,image.size(),betaBR_);
    regRotRG_ = cv::RotatedRect(center,image.size(),betaRG_);

    boundRotBG_ = regRotBG_.boundingRect2f();
    boundRotBR_ = regRotBR_.boundingRect2f();
    boundRotRG_ = regRotRG_.boundingRect2f();

    rotBG_ = cv::getRotationMatrix2D(center,betaBG_,1.0),
    rotBR_ = cv::getRotationMatrix2D(center,betaBR_,1.0),
    rotRG_ = cv::getRotationMatrix2D(center,betaRG_,1.0);
    rotBG_.at<double>(0,2) += 0.5*boundRotBG_.width  - 0.5*image.cols;
    rotBG_.at<double>(1,2) += 0.5*boundRotBG_.height - 0.5*image.rows;
    rotBR_.at<double>(0,2) += 0.5*boundRotBR_.width  - 0.5*image.cols;
    rotBR_.at<double>(1,2) += 0.5*boundRotBR_.height - 0.5*image.rows;
    rotRG_.at<double>(0,2) += 0.5*boundRotRG_.width  - 0.5*image.cols;
    rotRG_.at<double>(1,2) += 0.5*boundRotRG_.height - 0.5*image.rows;

  }

  // Record the image dimensions
  M_ = image.rows; N_ = image.cols; 

  end = clock();
  prepTime_ += ((double) end-begin)/CLOCKS_PER_SEC; 

}

/* ************************************************************************* */
void TriangleDisparity::setCameraCalib(double aBG, double aBR, double aRG,
  double cBG, double cBR, double cRG, double zDBG, double zDBR, double zDRG, 
  double sPx, double bBG, double bBR, double bRG, double fEff, double cGx, 
  double cGy, double gEff) { 
  

  aBG_ = aBG; aBR_ = aBR; aRG_ = aRG; 
  cBG_ = cBG; cBR_ = cBR; cRG_ = cRG; 
  zDBG_ = zDBG; zDBR_ = zDBR; zDRG_ = zDRG; 
  sPx_ = sPx; 
  bBG_ = bBG; bBR_ = bBR; bRG_ = bRG;
  fEff_ = fEff; cGx_ = cGx; cGy_ = cGy; gEff_ = gEff;
  isCamCalibSet_ = true;

}

/* ************************************************************************* */
void TriangleDisparity::calculateDisparity(cv::Mat& dispMapBG, 
  cv::Mat& dispMapBR, cv::Mat& dispMapRG, cv::Mat& validMap) { 

  clock_t begin = clock(), end;

  // Prepare all of the rotated images
  prepareImages();

  end = clock();
  prepTime_ += ((double) end-begin)/CLOCKS_PER_SEC; 
  begin = clock();

  // Create stereo block matcher and set minimum disparity
  // Use below for 1/4 scale 
  cv::Ptr<cv::StereoSGBM> sbm = cv::StereoSGBM::create(minDisparity_,
    numDisparities_,windowSize_,6*windowSize_*windowSize_,
    24*windowSize_*windowSize_,0,0,uniqueThresh_,800,1);

  /*
  cv::Ptr<cv::StereoBM> sbm = cv::StereoBM::create(numDisparities_,windowSize_);
  sbm->setMinDisparity(minDisparity_);
  sbm->setUniquenessRatio(uniqueThresh_);
  sbm->setTextureThreshold(0.0002); 
  */

  // Compute disparity (disparity is returned as a CV_16SC1 OpenCV matrix)
  // The value returned is 16x the actual disparity, and invalid values are 
  // one below the minimum disparity
  cv::Mat dispMapRotBG, dispMapRotBR, dispMapRotRG;
  sbm->compute(GrotBG_,BrotBG_,dispMapRotBG);
  sbm->compute(RrotBR_,BrotBR_,dispMapRotBR);
  sbm->compute(GrotRG_,RrotRG_,dispMapRotRG);

  // Rotate the disparity maps back to upright by reversing previous rotation
  // Use nearest neighbor so as not to average valid / invalid results
  cv::warpAffine(dispMapRotBG,dispMapBG,rotBG_,cv::Size(N_,M_),
                  cv::INTER_NEAREST + cv::WARP_INVERSE_MAP);
  cv::warpAffine(dispMapRotBR,dispMapBR,rotBR_,cv::Size(N_,M_),
                  cv::INTER_NEAREST + cv::WARP_INVERSE_MAP);
  cv::warpAffine(dispMapRotRG,dispMapRG,rotRG_,cv::Size(N_,M_),
                  cv::INTER_NEAREST + cv::WARP_INVERSE_MAP);

  // Convert to double and multiply by 1/16 to get disparity value
  dispMapBG.convertTo(dispMapBG,CV_64FC1,0.0625); 
  dispMapBR.convertTo(dispMapBR,CV_64FC1,0.0625);
  dispMapRG.convertTo(dispMapRG,CV_64FC1,0.0625);

  // Get shifted results for dispMapBR, moving them from red pinhole coords
  // (xR, yR) to green pinhole coords (xG, yG)
  double inv = minDisparity_ - 1.0;                     // invalid value
  cv::Mat dispMapBRs(M_,N_,CV_64FC1), 
          numAvgBRs(M_,N_,CV_8UC1);
  dispMapBRs = cv::Scalar(inv);
  numAvgBRs  = cv::Scalar(0);
  for (int yR = 0; yR < M_; yR++) { 

    double *BRptr  = dispMapBR.ptr<double>(yR);

    for (int xR = 0; xR < N_; xR++) { 

      // If disparity is valid, shift it from R coords to G coords by disparity
      if (BRptr[xR] != inv) { 

        int xG = xR + std::round(BRptr[xR]*pRtoG_.x),
            yG = yR + std::round(BRptr[xR]*pRtoG_.y);

        // Check that shifted coordinates are in the image; skip if not
        if (xG < 0 || xG > N_-1 || yG < 0 || yG > M_-1) { 
          continue;
        }

        char *numAvgPtr = &numAvgBRs.at<char>(yG,xG);
        (*numAvgPtr)++;

        dispMapBRs.at<double>(yG,xG) = (*numAvgPtr == 1)? BRptr[xR] : 
              dispMapBRs.at<double>(yG,xG) * ((double) *numAvgPtr-1.0)
                                           / ((double) *numAvgPtr    )
              + BRptr[xR]                  / ((double) *numAvgPtr    );
      }

    }
  }

  // Reassign dispMapBR to shifted version
  dispMapBR = dispMapBRs; 

  // Create validity map with binary flags in the least significant bits
  // Make lowest and highest disparity values invalid, since they are often
  // erroneous
  validMap.create(M_,N_,CV_8UC1);
  validMap = cv::Scalar(0);
  for (int y = 0; y < M_; y++) { 

    double *BGptr = dispMapBG.ptr<double>(y), *BRptr = dispMapBR.ptr<double>(y), 
           *RGptr = dispMapRG.ptr<double>(y);
    uchar  *Vptr  = validMap.ptr<uchar>(y);

    for (int x = 0; x < N_; x++) { 

      if (BGptr[x] <= minDisparity_ || BGptr[x] >= maxDisparity_) {
        BGptr[x] = inv;
      } 
      if (BRptr[x] <= minDisparity_ || BRptr[x] >= maxDisparity_) {
        BRptr[x] = inv;
      } 
      if (RGptr[x] <= minDisparity_ || RGptr[x] >= maxDisparity_) {
        RGptr[x] = inv; 
      }

      Vptr[x] = ( (BGptr[x] != inv)? 1 : 0 )    // bin: 0 0 0 0 0 0 0 (1 or 0)
              + ( (BRptr[x] != inv)? 2 : 0 )    // bin: 0 0 0 0 0 0 (1 or 0) 0
              + ( (RGptr[x] != inv)? 4 : 0 );   // bin: 0 0 0 0 0 (1 or 0) 0 0

    }

  }

  end = clock();
  calcTime_ += ((double) end-begin)/CLOCKS_PER_SEC; 

}

/* ************************************************************************* */
void TriangleDisparity::calculateDepth(const cv::Mat& dispMapBG, 
  const cv::Mat& dispMapBR, const cv::Mat& dispMapRG, const cv::Mat& validMap, 
  cv::Mat& depthMap) {

  clock_t begin = clock(), end;

  // Create a depth map with an invalid depth (-1.0 meter) by default
  depthMap.create(M_,N_,CV_64FC1);
  depthMap = cv::Scalar(-1.0);

  // Loop through all pixels
  for (int y = 0; y < M_; y++) { 

    const double 
          *BGptr = dispMapBG.ptr<double>(y), *BRptr = dispMapBR.ptr<double>(y), 
          *RGptr = dispMapRG.ptr<double>(y);
    const uchar *Vptr  = validMap.ptr<uchar>(y);
    double      *Dptr  = depthMap.ptr<double>(y);

    for (int x = 0; x < M_; x++) { 

      // Average valid results, or leave as invalid (-1.0) if no valid results
      double  depthBG     = zDBG_/(BGptr[x]*sPx_/bBG_ - 1.0 + aBG_)  - cBG_,
              depthBR     = zDBR_/(BRptr[x]*sPx_/bBR_ - 1.0 + aBR_)  - cBR_,
              depthRG     = zDRG_/(RGptr[x]*sPx_/bRG_ - 1.0 + aRG_)  - cRG_,
              BGvalid     = (double) ( Vptr[x] & 1      ),
              BRvalid     = (double) ((Vptr[x] & 2) >> 1),
              RGvalid     = (double) ((Vptr[x] & 4) >> 2),
              numValid    = BGvalid + BRvalid + RGvalid;
  
      if (numValid > 0.0) { 
        Dptr[x] = 0.0;
        Dptr[x] += (BGvalid == 1.0)? depthBG : 0.0;
        Dptr[x] += (BRvalid == 1.0)? depthBR : 0.0;
        Dptr[x] += (RGvalid == 1.0)? depthRG : 0.0;
        Dptr[x] /= numValid;
      }

    }
  }

  end = clock();
  calcTime_ += ((double) end-begin)/CLOCKS_PER_SEC; 

}

/* ************************************************************************* */
void TriangleDisparity::writeMap(const std::string& name, const cv::Mat& map,
                                    double minVal, double maxVal) { 

  std::string resDir("");

  // Get map dimensions
  int M = map.rows, N = map.cols;

  // Write the exact values of map to csv files
  std::string   mapFileName = resDir + name + ".csv";
  std::ofstream mapFile;
  mapFile.open(mapFileName.c_str());
  for (int i = 0; i < M; i++) { 
    const double *mapPtr = map.ptr<double>(i);
    for (int j = 0; j < N; j++) { 
      std::string comma = (j != N-1)? "," : ""; 
      mapFile << mapPtr[j] << comma;
    }
    mapFile << std::endl;
  }
  mapFile.close();

  // Convert map to type CV_8UC3 colormap with maximum contrast
  double    min, max, mu = (double) cv::mean(map)[0];
  cv::Mat   mapScaled8u, colorMap, validMask;
  cv::minMaxLoc(map, &min, &max);
  minVal = (minVal == -1000.0)? min : minVal;
  maxVal = (maxVal ==  1000.0)? max : maxVal;
  cv::subtract(map,cv::Scalar(minVal),map);
  map.convertTo(mapScaled8u, CV_8UC1, 1.0/(maxVal-minVal)*255.0);
  cv::applyColorMap(mapScaled8u,colorMap,cv::COLORMAP_JET);
  cv::imwrite(resDir + name + ".png",colorMap);

  // Print statistics
  std::cout << "Wrote " << resDir << name << ".csv" << std::endl << "  and " 
    << resDir << name << ".png" << std::endl;
  std::cout << "Min val: " << min<< ", Max val: " << max << 
    ", Mean val: " << mu << std::endl;

}

/* ************************************************************************* */
void TriangleDisparity::writePly(const std::string& name, 
  const cv::Mat& depthMap, double resizeFactor) {

  std::string resDir("");

  // Get map dimensions
  int M = depthMap.rows, N = depthMap.cols;
  if (M != M_ || N != N_) { 
    std::cout << "Error in TriangleDisparity::writePly Depth map image "
      << "dimensions must match input image dimensions" << std::endl;
    return;
  }

  // Determine number of valid pixels
  int numValid = 0;
  for (int i = 0; i < M; i++) { 
    const double *mapPtr = depthMap.ptr<double>(i);
    for (int j = 0; j < N; j++) { 
      if (mapPtr[j] != -1.0) { 
        numValid++;
      } 
    }
  }

  // Write the header
  std::string   mapFileName = resDir + name + ".ply";
  std::ofstream mapFile;
  mapFile.open(mapFileName.c_str());
  mapFile << "ply" << std::endl
          << "format ascii 1.0" << std::endl
          << "comment Generated by TriangleDisparity.cpp" << std::endl
          << "element vertex " << numValid << std::endl
          << "property float x" << std::endl
          << "property float y" << std::endl
          << "property float z" << std::endl
          << "property uchar red" << std::endl
          << "property uchar green" << std::endl
          << "property uchar blue" << std::endl
          << "end_header" << std::endl;

  // Write the 3D data, using the camera calibration, and using the green 
  // channel to provide grayscale color
  for (int i = 0; i < M; i++) { 
    const double *mapPtr = depthMap.ptr<double>(i);
    const uchar *Gptr = G_.ptr<uchar>(i);
    for (int j = 0; j < N; j++) { 
      if (mapPtr[j] != -1.0) { 
        double  z = mapPtr[j], 
                x = (j+1-cGx_*resizeFactor)*sPx_*(z-gEff_)/fEff_, 
                y = (i+1-cGy_*resizeFactor)*sPx_*(z-gEff_)/fEff_;
        int     I = std::min(255,3*(int) Gptr[j]);
        mapFile << x << " " << y << " " << z << " " << I << " " << I
                << " " << I << std::endl;
      }
    }
  }
  mapFile.close();

  std::cout << "Wrote " << resDir << name << ".ply" << std::endl;

}

/* ************************************************************************* */
void TriangleDisparity::writeImages(const std::string& name) { 
  
  std::string resDir("");

  // Convert to type CV_8U
  cv::imwrite(resDir + name + "-R.png",R_);
  cv::imwrite(resDir + name + "-G.png",G_);
  cv::imwrite(resDir + name + "-B.png",B_);
  cv::imwrite(resDir + name + "-BrotBG.png",BrotBG_);
  cv::imwrite(resDir + name + "-GrotBG.png",GrotBG_);
  cv::imwrite(resDir + name + "-BrotBR.png",BrotBR_);
  cv::imwrite(resDir + name + "-RrotBR.png",RrotBR_);
  cv::imwrite(resDir + name + "-RrotRG.png",RrotRG_);
  cv::imwrite(resDir + name + "-GrotRG.png",GrotRG_);

  std::cout << "Wrote 9 debugging images " << resDir << name << "-*.png" 
    << std::endl;

}

/* ************************************************************************* */
void TriangleDisparity::print() { 

  std::cout << "TriangleDisparity Properties ===================" << std::endl;
  std::cout << "  Transmissivity R   : " << tR_ << std::endl;
  std::cout << "  Transmissivity G   : " << tG_ << std::endl;
  std::cout << "  Transmissivity B   : " << tB_ << std::endl;
  std::cout << "  Pinhole centroid R : " << pRc_.x <<","<< pRc_.y << std::endl;
  std::cout << "  Pinhole centroid G : " << pGc_.x <<","<< pGc_.y << std::endl;
  std::cout << "  Pinhole centroid B : " << pBc_.x <<","<< pBc_.y << std::endl;
  std::cout << "  Crosstalk mat TCtoF: " << std::endl;
  std::cout << "    " << TCtoF_.at<double>(0,0) <<","<< TCtoF_.at<double>(0,1) 
            <<    "," << TCtoF_.at<double>(1,2) << std::endl;
  std::cout << "    " << TCtoF_.at<double>(1,0) <<","<< TCtoF_.at<double>(1,1) 
            <<    "," << TCtoF_.at<double>(0,2) << std::endl;
  std::cout << "    " << TCtoF_.at<double>(2,0) <<","<< TCtoF_.at<double>(2,1) 
            <<    "," << TCtoF_.at<double>(2,2) << std::endl;
  std::cout << "  Rotation angle BG  : " << betaBG_ << std::endl;
  std::cout << "  Rotation angle BR  : " << betaBR_ << std::endl;
  std::cout << "  Rotation angle RG  : " << betaRG_ << std::endl;
  std::cout << "  minDisparity       : " << minDisparity_ << std::endl;
  std::cout << "  maxDisparity       : " << maxDisparity_ << std::endl;
  std::cout << "  numDisparities     : " << numDisparities_ << std::endl;
  std::cout << "  uniqueThresh       : " << uniqueThresh_ << std::endl;
  std::cout << "  windowSize         : " << windowSize_ << std::endl;
  std::cout << "  Image height (M)   : " << M_ << std::endl;
  std::cout << "  Image width  (N)   : " << N_ << std::endl;
  std::cout << "Camera Calibration =============================" << std::endl;
  std::cout << "  aBG                : " << aBG_ << std::endl;
  std::cout << "  aBR                : " << aBR_ << std::endl;
  std::cout << "  aRG                : " << aRG_ << std::endl;
  std::cout << "  cBG                : " << cBG_ << std::endl;
  std::cout << "  cBR                : " << cBR_ << std::endl;
  std::cout << "  cRG                : " << cRG_ << std::endl;
  std::cout << "  zDBG               : " << zDBG_ << std::endl;
  std::cout << "  zDBR               : " << zDBR_ << std::endl;
  std::cout << "  zDRG               : " << zDRG_ << std::endl;
  std::cout << "  sPx                : " << sPx_ << std::endl;
  std::cout << "  bBG                : " << bBG_ << std::endl;
  std::cout << "  bBR                : " << bBR_ << std::endl;
  std::cout << "  bRG                : " << bRG_ << std::endl;
  std::cout << "  fEff               : " << fEff_ << std::endl;
  std::cout << "  cGx                : " << cGx_ << std::endl;
  std::cout << "  cGy                : " << cGy_ << std::endl;
  std::cout << "  gEff               : " << gEff_ << std::endl;
  std::cout << "================================================" << std::endl;

}

/* ************************************************************************* */
void TriangleDisparity::printTimes() { 

  std::cout << "TriangleDisparity Time --------" << std::endl;
  std::cout << "Preparation time [s] : " << prepTime_ << std::endl;
  std::cout << "Calculation time [s] : " << calcTime_ << std::endl;
  std::cout << "Total time [s]       : " << (prepTime_+calcTime_) 
                                                     << std::endl;
  std::cout << "-------------------------------" << std::endl;
}

/* ************************************************************************* */
void TriangleDisparity::resetTimes() { 

  prepTime_   = 0.0;  calcTime_   = 0.0;

}


/* ************************************************************************* */
void TriangleDisparity::prepareImages() { 

  // Rotate the images
  cv::warpAffine(B_,BrotBG_,rotBG_,boundRotBG_.size(),cv::INTER_LINEAR);
  cv::warpAffine(G_,GrotBG_,rotBG_,boundRotBG_.size(),cv::INTER_LINEAR);
  cv::warpAffine(B_,BrotBR_,rotBR_,boundRotBR_.size(),cv::INTER_LINEAR);
  cv::warpAffine(R_,RrotBR_,rotBR_,boundRotBR_.size(),cv::INTER_LINEAR);
  cv::warpAffine(R_,RrotRG_,rotRG_,boundRotRG_.size(),cv::INTER_LINEAR);
  cv::warpAffine(G_,GrotRG_,rotRG_,boundRotRG_.size(),cv::INTER_LINEAR);
 
}
