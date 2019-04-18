CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

PROJECT(CodedAperture)

# Choose between Debug or Release mode
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)

# For both Debug and Release modes
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -lrt -lpthread -pthread")

# For Debug mode
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -std=c++0x -g3 -O0")
# For Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -std=c++0x -O3")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funsafe-loop-optimizations -fsee -funroll-loops -fno-math-errno -funsafe-math-optimizations -ffinite-math-only -fno-signed-zeros") 

# Find packages
LINK_DIRECTORIES(/usr/local/lib)
INCLUDE_DIRECTORIES(/usr/local/include)
FIND_PACKAGE(OpenCV REQUIRED)

# Create triangle-disparity and triangle-disparity-static libraries
ADD_LIBRARY(triangle-disparity SHARED TriangleDisparity.cpp)
ADD_LIBRARY(triangle-disparity-static STATIC TriangleDisparity.cpp)
TARGET_LINK_LIBRARIES(triangle-disparity ${OpenCV_LIBS}) 
TARGET_LINK_LIBRARIES(triangle-disparity-static ${OpenCV_LIBS}) 
INSTALL(
  FILES TriangleDisparity.h 
  DESTINATION "${PROJECT_SOURCE_DIR}/include"
)
INSTALL(
  TARGETS triangle-disparity triangle-disparity-static
  LIBRARY DESTINATION "${PROJECT_SOURCE_DIR}/lib"
  ARCHIVE DESTINATION "${PROJECT_SOURCE_DIR}/lib"
) 

# Add tests directory
ADD_SUBDIRECTORY(tests)
