# Allow headers in tests to be included like
# #include "Histogram.h" rather than needing
# #include "Histogram/Histogram.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestHistogram TestHistogram.cpp)
target_link_libraries(TestHistogram)
