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

// STL
#include <numeric> // for 'accumulate'

// ITK
#include "itkNthElementImageAdaptor.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"

// Submodules
#include <ITKHelpers/ITKHelpers.h>

template <typename TBinValue, typename TQuadrantProperties>
template <typename TImage>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType HistogramGenerator<TBinValue, TQuadrantProperties>::ComputeScalarImageHistogram(
    const TImage* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
    const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax,
    const bool allowOutside)
{
  // Compute the histogram for each channel separately
  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    std::vector<typename TImage::PixelType> pixelValues =
        ITKHelpers::GetPixelValuesInRegion(image, region);

    HistogramType histogram;
    if(allowOutside)
    {
      histogram = ScalarHistogramAllowOutside(pixelValues, numberOfBinsPerDimensions, rangeMin, rangeMax);
    }
    else
    {
      histogram = ScalarHistogram(pixelValues, numberOfBinsPerDimensions, rangeMin, rangeMax);
    }

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue, typename TQuadrantProperties>
template <typename TScalarImage>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType HistogramGenerator<TBinValue, TQuadrantProperties>::ComputeImageHistogram1D(
    const TScalarImage* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMin,
    const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMax,
    const bool allowOutside)
{
  return ComputeScalarImageHistogram(image, region, numberOfBinsPerDimensions, rangeMin, rangeMax, allowOutside);
}

template <typename TBinValue, typename TQuadrantProperties>
template <typename TComponent, unsigned int Dimension>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType HistogramGenerator<TBinValue, TQuadrantProperties>::ComputeImageHistogram1D(
    const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const TComponent& rangeMin,
    const TComponent& rangeMax,
    const bool allowOutside)
{
  assert(image);
  //return Compute1DConcatenatedHistogramOfMultiChannelImage(image, region, numberOfBinsPerDimensions, rangeMin, rangeMax);

  // For Image<CovariantVector>, we must use a NthElementImageAdaptor

  typedef itk::Image<itk::CovariantVector<TComponent, Dimension>, 2> ImageType;

  // Compute the histogram for each channel separately
  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel
    typedef itk::Image<typename TypeTraits<typename ImageType::PixelType>::ComponentType, 2> ScalarImageType;

    typedef itk::NthElementImageAdaptor<ImageType, typename ScalarImageType::PixelType> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SelectNthElement(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    HistogramType histogram = ComputeScalarImageHistogram(adaptor.GetPointer(), region, numberOfBinsPerDimensions, rangeMin, rangeMax, allowOutside);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue, typename TQuadrantProperties>
template <typename TComponent>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType HistogramGenerator<TBinValue, TQuadrantProperties>::ComputeImageHistogram1D(
    const itk::VectorImage<TComponent, 2>* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const TComponent& rangeMin,
    const TComponent& rangeMax,
    const bool allowOutside)
{
  // For VectorImage, we must use a VectorImageToImageAdapter
  //return Compute1DConcatenatedHistogramOfMultiChannelImage(image, region, numberOfBinsPerDimensions, rangeMin, rangeMax);

  typedef itk::VectorImage<TComponent, 2> ImageType;
  // Compute the histogram for each channel separately
  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel
    typedef itk::Image<typename TypeTraits<typename ImageType::PixelType>::ComponentType, 2> ScalarImageType;

    typedef itk::VectorImageToImageAdaptor<typename TypeTraits<typename ImageType::PixelType>::ComponentType, 2> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SetExtractComponentIndex(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    HistogramType histogram = ComputeScalarImageHistogram(adaptor.GetPointer(), region, numberOfBinsPerDimensions, rangeMin, rangeMax, allowOutside);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue, typename TQuadrantProperties>
template <typename TValue>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType
HistogramGenerator<TBinValue, TQuadrantProperties>::ScalarHistogramAllowOutside(const std::vector<TValue>& values,
                                                                                const unsigned int numberOfBins,
                                                                                const TValue& rangeMin, const TValue& rangeMax)
{
  assert(numberOfBins > 0);

  assert(rangeMax > rangeMin);

  // Count how many values fall in each bin. We store these counts as floats because sometimes we want to normalize the counts.
  // std::cout << "Create histogram with " << numberOfBins << " bins." << std::endl;
  HistogramType bins(numberOfBins, 0);

  const float binWidth = (rangeMax - rangeMin) / static_cast<float>(numberOfBins);

  // If the bins are not reasonably sized, return an all zero histogram
  if(fabs(binWidth - 0.0f) < 1e-6)
  {
    return bins;
  }

  for(unsigned int i = 0; i < values.size(); ++i)
  {
    int bin = (values[i] - rangeMin) / binWidth;
    if(bin < 0)
    {
      bin = 0;
    }
    else if(bin >= static_cast<int>(numberOfBins)) // There are only (numberOfBins - 1) indexes since the bin ids start at 0
    {
      bin = numberOfBins - 1;
    }

    bins[bin]++;
  }

  return bins;
}

template <typename TBinValue, typename TQuadrantProperties>
template <typename TValue>
typename HistogramGenerator<TBinValue, TQuadrantProperties>::HistogramType HistogramGenerator<TBinValue, TQuadrantProperties>::ScalarHistogram(const std::vector<TValue>& values,
                                                                                                     const unsigned int numberOfBins,
                                                                                                     const TValue& rangeMin, const TValue& rangeMax)
{
  assert(numberOfBins > 0);

  assert(rangeMax > rangeMin);

  // Count how many values fall in each bin. We store these counts as floats because sometimes we want to normalize the counts.
  // std::cout << "Create histogram with " << numberOfBins << " bins." << std::endl;
  HistogramType bins(numberOfBins, 0);

  const float binWidth = (rangeMax - rangeMin) / static_cast<float>(numberOfBins);

  // If the bins are not reasonably sized, return an all zero histogram
  if(fabs(binWidth - 0.0f) < 1e-6)
  {
    return bins;
  }

  for(unsigned int i = 0; i < values.size(); ++i)
  {
    // Handle the special case of the value being exactly equal to the top of the range.
    // We cannot simply add epsilon fudge factors because in a standard case of uchar values (pixels),
    // 255 is a valid value, and 255+epsilon overflows a uchar.
    if(values[i] == rangeMax)
    {
      bins[bins.size() - 1]++; // Add 1 to the last bin
      continue;
    }

    int bin = (values[i] - rangeMin) / binWidth;
    if(bin < 0)
    {
      std::stringstream ss;
      ss << "Can't write to bin " << bin << "!" << std::endl;
      ss << "There are " << values.size() << " values." << std::endl;
      ss << "Range min " << static_cast<float>(rangeMin) << std::endl;
      ss << "Range max " << static_cast<float>(rangeMax) << std::endl;
      ss << "values[i] (i = " << i << ") = " << static_cast<float>(values[i]) << std::endl;
      ss << "binWidth " << binWidth << std::endl;
      throw std::runtime_error(ss.str());
    }
    else if(bin >= static_cast<int>(numberOfBins)) // There are only (numberOfBins - 1) indexes since the bin ids start at 0
    {
      std::stringstream ss;
      ss << "Can't write to bin " << bin << "!" << std::endl;
      ss << "There are " << values.size() << " values." << std::endl;
      ss << "Range min " << static_cast<float>(rangeMin) << std::endl;
      ss << "Range max " << static_cast<float>(rangeMax) << std::endl;
      ss << "values[i] (i = " << i << ") = " << static_cast<float>(values[i]) << std::endl;
      ss << "binWidth " << binWidth << std::endl;
      throw std::runtime_error(ss.str());
    }
    else // all is ok
    {
      bins[bin]++;
    }
  }

  return bins;
}

/*
template <typename TBinValue>
typename HistogramGenerator<TBinValue>::HistogramType HistogramGenerator<TBinValue>::ComputeHistogramOfGradient(const FloatVector2ImageType* gradientImage, const itk::ImageRegion<2>& region)
{
  // Discretize the continuum of possible angles of vectors in the right half-plane.
  // (We flip all vectors so that they are pointing towards the right half-plane, since their orientation is not important).
  // Add the strength/magnitude of the vector to the appropriate bin in the histogram.

  // The synatx is atan2(y,x)
  //float radiansToDegrees = 180./3.14159;
  float minValue = atan2(1, 0);
  float maxValue = atan2(-1, 0);

  float range = maxValue - minValue;

  unsigned int numberOfBins = 20;

  HistogramType histogram(numberOfBins);

  float step = range / static_cast<float>(numberOfBins);

  itk::ImageRegionConstIterator<FloatVector2ImageType> iterator(gradientImage, region);

  while(!iterator.IsAtEnd())
    {
    FloatVector2ImageType::PixelType v = iterator.Get();
    if(v[0] < 0)
      {
      v *= -1.0f;
      }
    float angle = atan2(v[1], v[0]);
    unsigned int bin = static_cast<unsigned int>( (angle - minValue) / step );
    histogram[bin] += iterator.Get().GetNorm();
    ++iterator;
    }

  // Normalize
  Helpers::NormalizeVector(histogram);

  return histogram;
}*/


#endif
