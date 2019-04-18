/**
 * @file testTriangleDisparity
 * @brief Test the TriangleDisparity class on a color-coded aperture image
 * @author Timothy P Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)
 */

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <math.h>
#include <stdlib.h>
#include "opencv2/opencv.hpp"
#include "opencv2/highgui.hpp"
#include "../TriangleDisparity.h"

using namespace cv;

const char *windowDisparity = "Disparity";

/**
 * @function main
 * @brief Main function
 */
int main(int argc, char** argv)
{

  // Deal with case where not enough arguments are passed
  if (argc < 5) { 
    std::cout << "Process a defocused image from a color-coded aperture"
      << std::endl 
      << "Usage: " << std::endl
      << "  testTriangleDisparity apertureName imageName "
      <<                      "minDepth maxDepth" << std::endl
      << "    apertureName  : either 'glass50' or 'glass85'" << std::endl
      << "    imageName     : path to the image to process"  << std::endl
      << "    minDepth      : minimum depth to consider [m]" << std::endl
      << "    maxDepth      : maximum depth to consider [m]" << std::endl;
      return 1;   
  }

  std::string apertureName = argv[1], imageName = argv[2];
  double  minDepth = std::max(0.5,atof(argv[3])), maxDepth = atof(argv[4]),
          tR, tG, tB, aBG, aBR, aRG, cBG, cBR, cRG, zDBG, zDBR, zDRG, 
          bBG, bBR, bRG, fEff, cGx, cGy, gEff, sPx;
  cv::Mat TCtoF;
  cv::Point2f pRc, pGc, pBc;

  if (apertureName == "glass50") { 

    // Camera calibration calib-glass50-led-08.24.2018
    tR = 0.368918; tG = 0.244706; tB = 0.386376;
    pRc = cv::Point2f(0.000000,0.000000);
    pGc = cv::Point2f(-0.321091,0.563854);
    pBc = cv::Point2f(-0.725798,0.029852);
    aBG = 1.001148; aBR = 1.000369; aRG = 1.001049,
    cBG = -0.063799; cBR = -0.064970; cRG = -0.060415;
    zDBG = 0.003363; zDBR = 0.001110; zDRG = 0.003052; sPx = 4.500000e-06;
    bBG = 0.293332; bBR = 0.893171; bRG = 0.323857; 
    fEff = 0.054958; cGx = 2456.156591; cGy = 2409.236163; gEff = -0.020838;

    TCtoF = (Mat_<double>(3,3) << 
          1.259430503189212, -0.296170305552113, -0.042522437167015 ,
         -0.506915825715605,  1.402281343425590, -0.305889481616210 ,
         -0.083531675722062, -0.523441191930818,  1.183625763845906 );

  } else {

    // Camera calibration calib-glass85-led-08.24.2018
    tR = 0.358558; tG = 0.289431; tB = 0.352011;
    pRc = cv::Point2f(0.000000,0.000000); 
    pGc = cv::Point2f(-0.760614,0.648781);
    pBc = cv::Point2f(-0.979609,-0.310073);
    aBG = 1.000389; aBR = 1.000205; aRG = 1.000229;
    cBG = -0.048268; cBR = -0.059000; cRG = -0.059806;
    zDBG = 0.003177; zDBR = 0.001761; zDRG = 0.001910; sPx = 4.500000e-06;
    bBG = 0.735652; bBR = 1.310770; bRG = 1.207753; 
    fEff = 0.095356; cGx = 1863.453699; cGy = 3102.308784; gEff = -0.210120;

    TCtoF = (Mat_<double>(3,3) << 
          1.257036294580806, -0.302448707745906, -0.035049432869650 ,
          -0.532219148021015, 1.422040705546023, -0.339437934835568 ,
          -0.048405772799910, -0.529890758982390, 1.200347639643513 );

  }

  // Time the execution
  clock_t begin = clock();  

  // Set a resize factor for the image
  double resizeFactor = 1.0/4.0;
  sPx = sPx / resizeFactor;

  // Set the disparity search space, SAD window size, and uniqueness threshold
  double dispMaxBG = bBG/sPx * ( 1.0 - aBG + zDBG/(minDepth + cBG) ),
         dispMaxBR = bBR/sPx * ( 1.0 - aBR + zDBR/(minDepth + cBR) ),
         dispMaxRG = bRG/sPx * ( 1.0 - aRG + zDRG/(minDepth + cRG) ),
         dispMinBG = bBG/sPx * ( 1.0 - aBG + zDBG/(maxDepth + cBG) ),
         dispMinBR = bBR/sPx * ( 1.0 - aBR + zDBR/(maxDepth + cBR) ),
         dispMinRG = bRG/sPx * ( 1.0 - aRG + zDRG/(maxDepth + cRG) ),
         dispMax = std::max( dispMaxBG, std::max( dispMaxBR, dispMaxRG ) ),
         dispMin = std::max( dispMinBG, std::max( dispMinBR, dispMinRG ) );
  int dispMax16       = (int) std::ceil( dispMax / 16.0 ), 
      dispMin16       = (int) std::floor( dispMin / 16.0 ),
      minDisparity    = dispMin16 * 16, 
      numDisparities  = dispMax16 * 16 - dispMin16 * 16,
      windowSize      = 15,   // 15 for 1/4 scale
      uniqueThresh    = 25;
  if (windowSize%2 == 0) {    // window size must be odd - fix if even
    windowSize++;
  }

  // Read the coded aperture image 
  Mat img8u, img8us;
  img8u = imread(argv[2]);

  // Resize the image, if necessary
  if (resizeFactor != 1.0) { 
    Size imDim(std::floor(resizeFactor*img8u.cols),
              std::floor(resizeFactor*img8u.rows));
    cv::resize(img8u,img8us,imDim,0,0,cv::INTER_LINEAR);
    std::cout << "Image resized to: " << img8us.rows << " x " << img8us.cols 
              << " pixels" << std::endl;
  } else {
    img8us = img8u;
    std::cout << "No image resize applied" << std::endl;
  }

  // Create the triangle disparity object
  TriangleDisparity triangleDisparity(tR,tG,tB,TCtoF,pRc,pGc,pBc,minDisparity,
    numDisparities,uniqueThresh,windowSize);

  // Set the camera calbration
  triangleDisparity.setCameraCalib(aBG,aBR,aRG,cBG,cBR,cRG,zDBG,zDBR,zDRG,
                                    sPx,bBG,bBR,bRG,fEff,cGx,cGy,gEff);

  // Disparity maps, validity map, and depth map
  Mat dispMapBG, dispMapBR, dispMapRG, validMap, depthMap;

  // Reset the computation times
  triangleDisparity.resetTimes();

  // Set the image, calculate the disparity, and write the results
  triangleDisparity.setImage(img8us);
  triangleDisparity.print();
  triangleDisparity.calculateDisparity(dispMapBG,dispMapBR,dispMapRG,validMap);

  // Calculate the depth
  triangleDisparity.calculateDepth(dispMapBG,dispMapBR,dispMapRG,
                                    validMap, depthMap);

  // Print times
  triangleDisparity.printTimes();

  // Write the results
  //triangleDisparity.writeMap(imageName + "-dispMapBG",dispMapBG);
  //triangleDisparity.writeMap(imageName + "-dispMapBR",dispMapBR);
  //triangleDisparity.writeMap(imageName + "-dispMapRG",dispMapRG);
  triangleDisparity.writeMap(imageName+"-depthMap",depthMap,minDepth,maxDepth);
  //triangleDisparity.writePly(imageName,depthMap,resizeFactor);

  // Write the debugging images
  triangleDisparity.writeImages(imageName);

  // Get execution time
  clock_t end = clock();
  double duration = ((double) end-begin)/CLOCKS_PER_SEC; 
  std::cout << "-------------------------------" << std::endl;
  std::cout << "testTriangleDisparity execution time: "
    << duration << "s" << std::endl;

  return 0;
}

