# MATLAB Software

This software implements algorithms described in:  
[1] T.P. Setterfield, A. Teran, N. Sigrist, D.B. Natzic, C.C. Liebe "Depth from
a Calibrated Color-Coded Aperture Camera using Disparity", 2019  
If using this software in a publication, please cite this paper.

## Pre-requisites

The third party PROPER library is required to generate point spread functions
(PSFs). The latest version can be found on the 
[PROPER SourceForge page](https://sourceforge.net/projects/proper-library/files/).
Be sure to download the MATLAB release, not the Python release.

Download this library and extract the contents to `thirdparty/proper`.

## Testing

To test generating color-coded aperture images at a fixed distance and
calculating their dense disparity, open MATLAB and run
`testTriangleDisparityProper.m`. Five different example test images from
the paper are available to choose from.

To test calculating depth from a real color-coded aperture image, open MATLAB
and run `testTriangleDisparityReal.m`.

To test the depth vs. disparity calibration procedure, open MATLAB and run
`calibrateTriangleDisparity.m`. To reduce the size of this repository, the raw
calibration data has been excluded. However, the previously processed data can
be loaded by setting `loadData = 1`.

To test the color space calibration procedure, open MATLAB and run 
`calibrateTriangleColorspace.m`.
