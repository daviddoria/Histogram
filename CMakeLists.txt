cmake_minimum_required(VERSION 2.8)

project(Histogram)

# This is needed for shared_ptr and the trick using enable_if and
# if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
  include_directories(${ITKHelpers_includes})
endif()

if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(Histogram_includes ${ITKHelpers_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  message("Histogram submodule used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# ITK
FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG
             ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology)
INCLUDE(${ITK_USE_FILE})

#add_library(Histogram)
add_custom_target(Histogram Histogram.h)

SET(Histogram_BuildTests OFF CACHE BOOL "Build Histogram tests?")
if(Histogram_BuildTests)
  add_subdirectory(Tests)
endif(Histogram_BuildTests)
