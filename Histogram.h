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

// Submodules
#include "Helpers/TypeTraits.h"

namespace Histogram
{

typedef float BinValueType;

typedef std::vector<BinValueType> HistogramType;

/** Compute the histograms of each channel of an image, and concatentate them together to form
  * of a 1D histogram. */

template <typename TImage>
HistogramType Compute1DConcatenatedHistogramOfMultiChannelImage(
                const TImage* image, const itk::ImageRegion<2>& region,
                const unsigned int numberOfBinsPerDimension,
                const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
                const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax);

/** Compute the histogram of a collection of values. */
template <typename TValue>
HistogramType ScalarHistogram(const std::vector<TValue>& values, const unsigned int numberOfBins,
                              const TValue& rangeMin, const TValue& rangeMax);

float HistogramIntersection(const HistogramType& histogram1, const HistogramType& histogram2);

float HistogramDifference(const HistogramType& histogram1, const HistogramType& histogram2);

void WriteHistogram(const HistogramType& histogram, const std::string& filename);

void OutputHistogram(const HistogramType& histogram);

} // end namespace

#include "Histogram.hpp"

#endif
