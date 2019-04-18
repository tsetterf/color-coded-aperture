# C++ Software

This software implements algorithms described in:  
[1] T.P. Setterfield, A. Teran, N. Sigrist, D.B. Natzic, C.C. Liebe "Depth from
a Calibrated Color-Coded Aperture Camera using Disparity", 2019  
If using this software in a publication, please cite this paper.

## Pre-requisites

CMake and OpenCV are required. The software has been tested with with CMake
version 3.5.1 and OpenCV version 2.4.9.1, but should work with other versions
as well. These libraries can be installed using the following command.
```
sudo apt-get install cmake libopencv-dev
```

## Compiling and Installing

To compile, navigate to the root directory and run the following commands.
```
mkdir build
cd build
cmake ..
make
```

If the build is successful, install the library (to lib folder) and binary (to
bin folder) using the following command.
```
make install
```

## Testing

The binary `testTriangleDisparity` can be used to process single images. To 
try out the algorithm on an example image from [1], navigate to the project
root directory and run the following command.
```
./bin/testTriangleDisparity glass50 ../data/blocks-glass50-09.26.2018/Dalsa-18.jpg 1.5 3.0
```

The following files will be written in the data directory containing 
Dalsa-18.jpg:
- `Dalsa-18.jpg-depthMap.csv`: Dense depth map as a comma separated value file
  of doubles.
- `Dalsa-18.jpg-depthMap.png`: Dense depth map as a PNG image colored by depth.
- `Dalsa-18.jpg-R.png`: Input image red channel as a grayscale PNG image.
- `Dalsa-18.jpg-G.png`: Input image green channel as a grayscale PNG image.
- `Dalsa-18.jpg-B.png`: Input image blue channel as a grayscale PNG image.
- `Dalsa-18.jpg-RrotRG.png`: Input image red channel rotated for horizontal
  block matching between the red and green channels.
- `Dalsa-18.jpg-GrotRG.png`: Input image green channel rotated for horizontal
  block matching between the red and green channels.
- `Dalsa-18.jpg-RrotBR.png`: Input image red channel rotated for horizontal
  block matching between the blue and red channels.
- `Dalsa-18.jpg-BrotBR.png`: Input image blue channel rotated for horizontal
  block matching between the blue and red channels.
- `Dalsa-18.jpg-GrotBG.png`: Input image green channel rotated for horizontal
  block matching between the blue and green channels.
- `Dalsa-18.jpg-BrotBG.png`: Input image blue channel rotated for horizontal
  block matching between the blue and green channels.





