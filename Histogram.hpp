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

#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "Histogram.h"

// Custom
#include "Helpers/Helpers.h"

// STL
#include <numeric> // for 'accumulate'

// ITK
#include "itkRegionOfInterestImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace Histogram
{

template <typename TImage>
std::vector<float> Compute1DConcatenatedHistogramOfMultiChannelImage(
                      const FloatVectorImageType* image,
                      const itk::ImageRegion<2>& region,
                      const unsigned int numberOfBinsPerDimensions,
                      const TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
                      const TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax)
{
  // Compute the histogram for each channel separately
  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel
    
    HistogramType histogram(channelHistograms[0]->GetSize(0));
    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}


template <typename TValue>
std::vector<float> ScalarHistogram(const std::vector<TValue>& values, const unsigned int numberOfBins,
                                   const TValue& rangeMin, const TValue& rangeMax)
{
  // Count how many values fall in each bin. We store these counts as floats because sometimes we want to normalize the counts.
  // std::cout << "Create histogram with " << numberOfBins << " bins." << std::endl;
  std::vector<float> bins(numberOfBins, 0);

  const float binWidth = (rangeMax - rangeMin) / static_cast<float>(numberOfBins);

  if(fabs(binWidth - 0.0f) < 1e-6)
  {
    return bins;
  }

  for(unsigned int i = 0; i < values.size(); ++i)
  {
    int bin = (values[i] - rangeMin) / binWidth;
    if(bin < 0)
    {
      std::stringstream ss;
      ss << "Can't write to bin " << bin << "!" << std::endl;
      ss << "There are " << values.size() << " values." << std::endl;
      ss << "Range min " << rangeMin << std::endl;
      ss << "Range max " << rangeMax << std::endl;
      ss << "values[i] (i = " << i << ") = " << values[i] << std::endl;
      ss << "binWidth " << binWidth << std::endl;
      throw std::runtime_error(ss.str());
    }
    bins[bin]++;
  }

  return bins;
}

float HistogramIntersection(const std::vector<float>& histogram1, const std::vector<float>& histogram2)
{
  if(histogram1.size() != histogram2.size())
    {
    std::cerr << "Histograms must be the same size!" << std::endl;
    return 0;
    }

  float totalIntersection = 0.0f;
  for(unsigned int bin = 0; bin < histogram1.size(); ++bin)
    {
    // The casts to float are necessary other wise the integer division always ends up = 0 !
    float frequency1 = static_cast<float>(histogram1[bin]);
    float frequency2 = static_cast<float>(histogram2[bin]);
    //std::cout << "frequency1: " << frequency1 << std::endl;
    //std::cout << "frequency2: " << frequency2 << std::endl;
    float intersection = std::min(frequency1, frequency2);
    //std::cout << "intersection: " << intersection << std::endl;
    totalIntersection += intersection;
    }

  float totalFrequency = std::accumulate(histogram1.begin(), histogram1.begin() + histogram1.size(), 0.0f);
  //std::cout << "totalFrequency: " << totalFrequency << std::endl;
  float normalizedIntersection = totalIntersection / totalFrequency;
  //std::cout << "normalizedIntersection: " << normalizedIntersection << std::endl;
  return normalizedIntersection;
}

void WriteHistogram(const std::vector<float>& histogram, const std::string& filename)
{
  std::ofstream fout(filename.c_str());
  for(unsigned int i = 0; i < histogram.size(); ++i)
    {
    fout << histogram[i] << " ";
    }

  fout.close();
}

void OutputHistogram(const HistogramType& histogram)
{
  for(unsigned int i = 0; i < histogram.size(); ++i)
    {
    std::cout << histogram[i] << " ";
    }
}

} // end namespace

#endif
