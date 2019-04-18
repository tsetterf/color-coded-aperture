/**
 * @file TriangleDisparity.h
 * @brief Finds a dense depth map from color-coded aperture image
 * @date May 18, 2018
 * @author Timothy Setterfield (Timothy.P.Setterfield@jpl.nasa.gov)
 */

#ifndef TRIANGLEDISPARITY_H_ 
#define TRIANGLEDISPARITY_H_

#include <vector>
#include "opencv2/opencv.hpp"
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"

/**
 *  The TriangleDisparity class analyzes images formed using the 
 *  three-pinhole color-coded aperture mask below:
 *            O  tR*R            R: red filter
 *                               G: green filter
 *                               B: blue filter
 *   tB*B  O      O  tG*G
 * 
 *  This mask features three pinholes, with red, green, and blue filters with
 *  different transmissivites:
 *      tR: transmissivity of red pinhole
 *      tG: transmissivity of green pinhole
 *      tB: transmissivity of blue pinhole
 *  
 *  The disparity between the blue and green images, between the blue
 *  and red images, and between the red and green images is calculated. To 
 *  use conventional block matching algorithms, the image is rotated so that 
 *  the disparity is in the horizontal direction for each "left"-"right"
 *  comparison.
 */

class TriangleDisparity {

private:

  /** Transmissivity of each pinhole */
  double tR_, tG_, tB_;

  /** Crosstalk from filter to channel */
  cv::Mat TCtoF_;

  /** Pinhole centroid locations in units of aperture radius */
  cv::Point2f pRc_, pGc_, pBc_;

  /** Minimum and maximum disparity [px] */
  int minDisparity_, maxDisparity_;

  /** Number of disparities for block matching [px] (must be multiple of 16) */
  int numDisparities_;

  /** 
   *  Threshold of uniquness value 100*(SAD2/SADopt - 1), where SADopt is the 
   *  minimum sum of asbolute differences, and SAD2 is the second-best
   *  non-adjacent value [%].
   */
  int uniqueThresh_;

  /** SAD window size for block matcher [px] (must be odd) */
  int windowSize_;

  /** Whether the camera calibration has been set */
  bool isCamCalibSet_;

  /** 
   *  Camera calibration parameters required for calculating depth using model:
   *    disparity   = b (1 - a + zD/(depth + c))              [m]
   *    depth       = zD / (disparityPx*sPx/b - 1 + a) - c    [m]
   */
  double aBG_, aBR_, aRG_, cBG_, cBR_, cRG_, zDBG_, zDBR_, zDRG_, sPx_,
         bBG_, bBR_, bRG_; 

  /** 
   *  Camera calibration parameters required for projecting 3D points onto the 
   *  image plane using model (cGx, cGy are MATLAB-generated; base index is 1):  
   *    u = fEff*x/(z - gEff) + (cGx-1)                        [px]
   *    v = fEff*y/(z - gEff) + (cGy-1)                        [px]
   */
  double fEff_, cGx_, cGy_, gEff_;

  /** Blurred red, green, and blue image channels. */
  cv::Mat R_, G_, B_;

  /** Unit vectors between aperture centers in image coordinates */
  cv::Point2f pBtoG_, pBtoR_, pRtoG_;

  /** Rotation angles to facilitate horizontal disparity finding [deg] */
  double betaBG_, betaBR_, betaRG_;

  /** Affine transformation matrix for 2D rotations */
  cv::Mat rotBG_, rotBR_, rotRG_;

  /** Rotated versions of the red, green, and blue channels */
  cv::Mat BrotBG_, GrotBG_, BrotBR_, RrotBR_, RrotRG_, GrotRG_;

  /** Rotated rectangle encompassing rotated image */
  cv::RotatedRect regRotBG_, regRotBR_, regRotRG_;

  /** Rectangle forming outer bounding box around rotated image */
  cv::Rect2f boundRotBG_, boundRotBR_, boundRotRG_;

  /** Image dimensions M = # of rows, N = # of columns */
  int M_, N_; 

  /** Execution times for preparing images, and calculating disparity */
  double prepTime_, calcTime_;

public:

  /** Constructor */ 
  TriangleDisparity(double tR, double tG, double tB, cv::Mat TCtoF, 
    cv::Point2f pRc, cv::Point2f pGc, cv::Point2f pBc, int minDisparity, 
    int numDisparities, int uniqueThresh, int windowSize);

  /** Destructor */
  virtual ~TriangleDisparity();
  
  /** Set the input full image and extract red, green, blue channels. */
  void setImage(const cv::Mat& image);

  /** Set the camera calibration parameters */
  void setCameraCalib(double aBG, double aBR, double aRG, double cBG, 
    double cBR, double cRG, double zDBG, double zDBR, double zDRG, double sPx, 
    double bBG, double bBR, double bRG, double fEff, double cGx, double cGy,
    double gEff);

  /** 
   *  Calculate the disparity maps and return them with double precision 
   *  (CV_64FC1). Also return the validity map, which is a CV_8UC1 matrix
   *  containing a a binary flag in the three least significant bits, 
   *    i.e.: 0 0 0 0 0 isDispMapRGvalid isDispMapBRvalid isDispMapBGvalid 
   *    e.g.: 0 0 0 0 0 1 1 0 = 6, means dispMapRG and dispMapBG are valid.
   */
  void calculateDisparity(cv::Mat& dispMapBG, cv::Mat& dispMapBR, 
          cv::Mat& dispMapRG, cv::Mat& validMap);

  /** 
   *  Calculate the consensus depth map by averaging the valid results from the 
   *  existing disparity maps.
   */
  void calculateDepth(const cv::Mat& dispMapBG, const cv::Mat& dispMapBR, 
    const cv::Mat& dispMapRG, const cv::Mat& validMap, cv::Mat& depthMap);

  /** 
   *  Write the resulting disparity and depth map to a file; for visualization
      use minVal and maxVal for scaling.
   */
  void writeMap(const std::string& name, const cv::Mat& map, 
                  double minVal=-1000.0, double maxVal=1000.0);

  /** Write the point cloud in .ply format  to a file */
  void writePly(const std::string& name, const cv::Mat& depthMap,
                  double resizeFactor);

  /** Write relevant images to files for debugging */
  void writeImages(const std::string& name);

  /** Print information about the triangle disparity object */
  void print();

  /** Print a report of the computation times */
  void printTimes();

  /** Reset the compuation times */
  void resetTimes();

private:

  /** Prepare the rotated images */
  void prepareImages();

};

#endif /* TRIANGLEDISPARITY_H_ */
