# Allow headers in tests to be included like
# #include "Histogram.h" rather than needing
# #include "Histogram/Histogram.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestHistogram TestHistogram.cpp)
target_link_libraries(TestHistogram ${Histogram_libraries})

add_executable(TestMaskedHistogram TestMaskedHistogram.cpp)
target_link_libraries(TestMaskedHistogram ${Histogram_libraries})

add_executable(TestHistogramTiming TestHistogramTiming.cpp)
target_link_libraries(TestHistogramTiming ${Histogram_libraries})

message("Histogram_libraries in Histogram/Tests is ${Histogram_libraries}")
