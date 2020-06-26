cmake_minimum_required (VERSION 3.9)
project (TestC4Detector)

option(BUILD_SHARED_LIBS "Build shared libs" ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

#dont use find_package(OpenCV)
#(OpenCVConfig.cmake has some problems with conan package)
set(OpenCV_LIBS CONAN_PKG::opencv)

ADD_EXECUTABLE (TestC4Detector c4-pedestrian-detector.cpp)

find_package(OpenCV REQUIRED)

if (OpenCV_FOUND)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
#	include_directories("${OPENCV_CONFIG_PATH}/include")
	target_link_libraries( TestC4Detector ${OpenCV_LIBS} )
else (OpenCV_FOUND)
message("PB->OPENCV = ${OpenCV_INCLUDE_DIRS}")
endif (OpenCV_FOUND)
