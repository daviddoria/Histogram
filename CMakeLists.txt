cmake_minimum_required(VERSION 2.8)

project(Histogram)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# This is needed for shared_ptr and the trick using enable_if and
# if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology)
  INCLUDE(${ITK_USE_FILE})
endif()

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
  include_directories(${ITKHelpers_includes})
  set(Histogram_includes ${Histogram_includes} ${ITKHelpers_includes})
endif()

# Whether or not the submodule is added here, this library is required.
set(Histogram_libraries ${Histogram_libraries} ITKHelpers)

if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(Histogram_includes ${Histogram_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(Histogram_libraries ${Histogram_libraries} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  message("Histogram submodule used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

add_custom_target(Histogram SOURCES Histogram.h Histogram.hpp MaskedHistogram.h MaskedHistogram.hpp)

SET(Histogram_BuildTests OFF CACHE BOOL "Build Histogram tests?")
if(Histogram_BuildTests)
  add_subdirectory(Tests)
endif(Histogram_BuildTests)
