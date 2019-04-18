# Include the files in the source directory
INCLUDE_DIRECTORIES(../)

# TriangleDisparity
ADD_EXECUTABLE(testTriangleDisparity testTriangleDisparity.cpp)
TARGET_LINK_LIBRARIES(testTriangleDisparity ${OpenCV_LIBS} 
  triangle-disparity-static)
INSTALL(TARGETS testTriangleDisparity
  RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin"
) 

