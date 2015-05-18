cmake_minimum_required(VERSION 2.8)

project(Histogram)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

add_subdirectory(CMakeHelpers)

# Enable c++11
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology)
  INCLUDE(${ITK_USE_FILE})
endif()

UseSubmodule(ITKHelpers Histogram)

CreateSubmodule(Histogram)

# Add non-compiled files to the project
add_custom_target(Histogram SOURCES
Histogram.h
HistogramGenerator.h
HistogramGenerator.hpp
MaskedHistogramGenerator.h
MaskedHistogramGenerator.hpp
HistogramHelpers.hpp
HistogramDifferences.hpp
QuadrantHistogram.h)

# Build the tests if requested
option(Histogram_BuildTests "Build Histogram tests?" OFF)
if(Histogram_BuildTests)
  add_subdirectory(Tests)
endif(Histogram_BuildTests)

# Display where this code was used from (if it is used as a submodule, there may be multiple instances of this submodule in the project, only the first of which is used)
option(Histogram_ShowSubmoduleLocation "Show the path from which Histogram was used?" OFF)
if(Histogram_ShowSubmoduleLocation)
  message("Histogram used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif(Histogram_ShowSubmoduleLocation)
