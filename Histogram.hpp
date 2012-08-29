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

template <typename TBinValue>
template <typename TImage>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ComputeScalarImageHistogram(
    const TImage* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
    const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax)
{
  // Compute the histogram for each channel separately
  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    std::vector<typename TImage::PixelType> pixelValues =
        ITKHelpers::GetPixelValuesInRegion(image, region);

    HistogramType histogram = Histogram::ScalarHistogram(pixelValues, numberOfBinsPerDimensions, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue>
template <typename TScalarImage>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ComputeImageHistogram1D(
    const TScalarImage* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMin,
    const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMax)
{
  return ComputeScalarImageHistogram(image, region, numberOfBinsPerDimensions, rangeMin, rangeMax);
}

template <typename TBinValue>
template <typename TComponent, unsigned int Dimension>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ComputeImageHistogram1D(
    const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const TComponent& rangeMin,
    const TComponent& rangeMax)
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

    HistogramType histogram = ComputeScalarImageHistogram(adaptor.GetPointer(), region, numberOfBinsPerDimensions, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue>
template <typename TComponent>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ComputeImageHistogram1D(
    const itk::VectorImage<TComponent, 2>* image,
    const itk::ImageRegion<2>& region,
    const unsigned int numberOfBinsPerDimensions,
    const TComponent& rangeMin,
    const TComponent& rangeMax)
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

    HistogramType histogram = ComputeScalarImageHistogram(adaptor.GetPointer(), region, numberOfBinsPerDimensions, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue>
template <typename TValue>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ScalarHistogram(const std::vector<TValue>& values, const unsigned int numberOfBins,
                                                                                   const TValue& rangeMin, const TValue& rangeMax)
{
  // Count how many values fall in each bin. We store these counts as floats because sometimes we want to normalize the counts.
  // std::cout << "Create histogram with " << numberOfBins << " bins." << std::endl;
  Histogram::HistogramType bins(numberOfBins, 0);

  const float binWidth = (rangeMax - rangeMin) / static_cast<float>(numberOfBins);

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

template <typename TBinValue>
float Histogram<TBinValue>::HistogramCoherence(const HistogramType& idealHistogram, const HistogramType& queryHistogram)
{
  if(idealHistogram.size() != queryHistogram.size())
  {
    std::stringstream ss;
    ss << "HistogramCoherence: Histograms must be the same size! idealHistogram is " << idealHistogram.size() << " while queryHistogram is " << queryHistogram.size();
    throw std::runtime_error(ss.str());
  }

  unsigned int numberOfNewBins = 0;
  for(unsigned int bin = 0; bin < idealHistogram.size(); ++bin)
  {
    if((idealHistogram[bin] == 0) && (queryHistogram[bin] > 0))
    {
      numberOfNewBins++;
    }
  }

  return numberOfNewBins;
}

template <typename TBinValue>
float Histogram<TBinValue>::WeightedHistogramDifference(const HistogramType& idealHistogram, const HistogramType& queryHistogram)
{
  // assert(TBinValue is a signed type)
  if(idealHistogram.size() != queryHistogram.size())
  {
    std::stringstream ss;
    ss << "Histograms must be the same size! idealHistogram is " << idealHistogram.size()
       << " while queryHistogram is " << queryHistogram.size();
    throw std::runtime_error(ss.str());
  }

  float difference = 0.0f;

  for(unsigned int bin = 0; bin < idealHistogram.size(); ++bin)
  {
    // Large values get a small weight
    //float weight = 1.0f/static_cast<float>(idealHistogram[bin]); // cannot do this, because some bins are 0
    float weight = 100.0f - static_cast<float>(idealHistogram[bin]); // this 100 is arbitrary, but if our patches have only a few hundred
        // pixels, and there are ~50 bins, then there should never be a bin with > 100 count.
    difference += weight * fabs(static_cast<float>(idealHistogram[bin]) - static_cast<float>(queryHistogram[bin]));
  }

  return difference;
}

template <typename TBinValue>
float Histogram<TBinValue>::HistogramDifference(const HistogramType& histogram1, const HistogramType& histogram2)
{
  // assert(TBinValue is a signed type)
  if(histogram1.size() != histogram2.size())
  {
    std::stringstream ss;
    ss << "Histograms must be the same size! histogram1 is " << histogram1.size() << " while histogram2 is " << histogram2.size();
    throw std::runtime_error(ss.str());
  }

  //float difference = 0.0f;
  TBinValue difference = 0;
  for(unsigned int bin = 0; bin < histogram1.size(); ++bin)
  {
    //difference += fabs(static_cast<float>(histogram1[bin]) - static_cast<float>(histogram2[bin]));
    difference += fabs(histogram1[bin] - histogram2[bin]);
  }

  return difference;
}

template <typename TBinValue>
float Histogram<TBinValue>::HistogramIntersection(const HistogramType& histogram1, const HistogramType& histogram2)
{
  if(histogram1.size() != histogram2.size())
  {
    std::cerr << "Histograms must be the same size!" << std::endl;
    return 0;
  }

  TBinValue totalIntersection = 0.0f;
  for(unsigned int bin = 0; bin < histogram1.size(); ++bin)
  {
    // The casts to float are necessary other wise the integer division always ends up = 0 !
    TBinValue binCount1 = histogram1[bin];
    TBinValue binCount2 = histogram2[bin];
    //std::cout << "frequency1: " << frequency1 << std::endl;
    //std::cout << "frequency2: " << frequency2 << std::endl;
    TBinValue intersection = std::min(binCount1, binCount2);
    //std::cout << "intersection: " << intersection << std::endl;
    totalIntersection += intersection;
  }

  // Why is this of histogram1 instead of histogram2? Are we assuming they have the same total bin sum?
  float totalFrequency = std::accumulate(histogram1.begin(), histogram1.begin() + histogram1.size(), 0.0f); // This 0.0f is the "initial value of the sum". It has to be a float, or the sum will be an integer sum.
  //std::cout << "totalFrequency: " << totalFrequency << std::endl;
  float normalizedIntersection = static_cast<float>(totalIntersection) / totalFrequency;
  //std::cout << "normalizedIntersection: " << normalizedIntersection << std::endl;
  return normalizedIntersection;
}


template <typename TBinValue>
void Histogram<TBinValue>::WriteHistogram(const Histogram::HistogramType& histogram, const std::string& filename)
{
  std::ofstream fout(filename.c_str());
  for(unsigned int i = 0; i < histogram.size(); ++i)
  {
    fout << histogram[i] << " ";
  }

  fout.close();
}

template <typename TBinValue>
void Histogram<TBinValue>::OutputHistogram(const HistogramType& histogram)
{
  for(unsigned int i = 0; i < histogram.size(); ++i)
  {
    std::cout << histogram[i] << " ";
  }
}
/*
template <typename TBinValue>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::ComputeHistogramOfGradient(const FloatVector2ImageType* gradientImage, const itk::ImageRegion<2>& region)
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
