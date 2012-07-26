/*=========================================================================
 *
 *  Copyright David Doria 2012 daviddoria@gmail.com
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

// STL
#include <string>
#include <vector>

// ITK
#include "itkImageRegion.h"

namespace Histogram
{

template <typename TImage>
std::vector<float> Compute1DHistogramOfMultiChannelImage(const TImage* image,
                                                         const itk::ImageRegion<2>& region,
                                                         const unsigned int numberOfBins);


float HistogramIntersection(const std::vector<float>& histogram1, const std::vector<float>& histogram2);

void WriteHistogram(const std::vector<float>& histogram1, const std::string& filename);

} // end namespace

#endif
