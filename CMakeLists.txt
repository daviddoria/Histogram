cmake_minimum_required(VERSION 2.8)

project(Histogram)

if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(Histogram_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endif()

# This is needed for shared_ptr and the trick using enable_if and
# if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

if(NOT TARGET Helpers)
  add_subdirectory(Helpers)
  include_directories(${Helpers_includes})
endif()

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG
             ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology)
INCLUDE(${ITK_USE_FILE})

#add_library(Histogram Histogram.cpp)
# add_library(Histogram)

SET(Histogram_BuildTests OFF CACHE BOOL "Build Histogram tests?")
if(Histogram_BuildTests)
  add_subdirectory(Tests)
endif(Histogram_BuildTests)
