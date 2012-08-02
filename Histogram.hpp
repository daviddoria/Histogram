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
  //return Compute1DConcatenatedHistogramOfMultiChannelImage(image, region, numberOfBinsPerDimensions, rangeMin, rangeMax);

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
/*
template <typename TBinValue>
std::vector<typename Histogram<TBinValue>::HistogramType>
Histogram<TBinValue>::ComputeHistogramsOfRegion(const FloatVectorImageType* image, const itk::ImageRegion<2>& region)
{

  std::vector<HistogramType::Pointer> histograms;

  typedef itk::RegionOfInterestImageFilter< FloatVectorImageType, FloatVectorImageType > RegionOfInterestImageFilterType;
  RegionOfInterestImageFilterType::Pointer regionOfInterestImageFilter = RegionOfInterestImageFilterType::New();
  regionOfInterestImageFilter->SetRegionOfInterest(region);
  regionOfInterestImageFilter->SetInput(image);
  regionOfInterestImageFilter->Update();

  // Compute the histogram of each channel
  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {
    typedef itk::VectorIndexSelectionCastImageFilter<FloatVectorImageType, FloatScalarImageType> IndexSelectionType;
    IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->GetOutput()->SetRequestedRegion(region);
    indexSelectionFilter->Update();

    const unsigned int MeasurementVectorSize = 1;
    const unsigned int binsPerDimension = 30;

    HistogramType::MeasurementVectorType lowerBound(MeasurementVectorSize);
    lowerBound.Fill(0);

    HistogramType::MeasurementVectorType upperBound(MeasurementVectorSize);
    upperBound.Fill(255) ;

    HistogramType::SizeType size(MeasurementVectorSize);
    size.Fill(binsPerDimension);

    ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();
    imageToHistogramFilter->SetInput(indexSelectionFilter->GetOutput());
    imageToHistogramFilter->SetHistogramBinMinimum(lowerBound);
    imageToHistogramFilter->SetHistogramBinMaximum(upperBound);
    imageToHistogramFilter->SetHistogramSize(size);
    imageToHistogramFilter->SetAutoMinimumMaximum(false);
    imageToHistogramFilter->Update();

    histograms.push_back(imageToHistogramFilter->GetOutput());
    }

  return histograms;
}*/
/*
template <typename TBinValue>
std::vector<HistogramType::Pointer> Histogram<TBinValue>::ComputeHistogramsOfMaskedRegion(const FloatVectorImageType* image, const itk::ImageRegion<2>& imageRegion, const Mask* mask, const itk::ImageRegion<2>& maskRegion, const unsigned int binsPerDimension)
{

  std::vector<HistogramType::Pointer> histograms;

  // Compute the histogram of each channel
  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {

    typedef itk::VectorIndexSelectionCastImageFilter<FloatVectorImageType, FloatScalarImageType> IndexSelectionType;
    IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->GetOutput()->SetRequestedRegion(imageRegion);
    indexSelectionFilter->Update();

    const unsigned int MeasurementVectorSize = 1;

    HistogramType::MeasurementVectorType lowerBound(MeasurementVectorSize);
    lowerBound.Fill(0);

    HistogramType::MeasurementVectorType upperBound(MeasurementVectorSize);
    upperBound.Fill(255) ;

    HistogramType::SizeType size(MeasurementVectorSize);
    size.Fill(binsPerDimension);


    itk::ImageRegionConstIterator<FloatScalarImageType> imageIterator(indexSelectionFilter->GetOutput(), imageRegion);
    itk::ImageRegionConstIterator<Mask> maskIterator(mask, maskRegion);

    HistogramType::Pointer histogram = HistogramType::New();

    histogram->SetMeasurementVectorSize(MeasurementVectorSize);

    histogram->Initialize(size, lowerBound, upperBound );

    while(!imageIterator.IsAtEnd())
      {
      if(mask->IsHole(maskIterator.GetIndex()))
        {
        ++imageIterator;
        ++maskIterator;
        continue;
        }
      float pixel = imageIterator.Get();
      HistogramType::MeasurementVectorType mv(MeasurementVectorSize);
      mv[0] = pixel;

      histogram->IncreaseFrequencyOfMeasurement(mv, 1);

      ++imageIterator;
      ++maskIterator;
      }

    histograms.push_back(histogram);
    }

  return histograms;
}*/
/*
template <typename TBinValue>
std::vector<HistogramType::Pointer> Histogram<TBinValue>::ComputeHistogramsOfRegionManual(const FloatVectorImageType* image, const itk::ImageRegion<2>& region, const unsigned int numberOfBins)
{

  std::vector<HistogramType::Pointer> histograms;

  // Compute the histogram of each channel
  for(unsigned int i = 0; i < image->GetNumberOfComponentsPerPixel(); ++i)
    {

    typedef itk::VectorIndexSelectionCastImageFilter<FloatVectorImageType, FloatScalarImageType> IndexSelectionType;
    IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
    indexSelectionFilter->SetIndex(i);
    indexSelectionFilter->SetInput(image);
    indexSelectionFilter->GetOutput()->SetRequestedRegion(region);
    indexSelectionFilter->Update();

    const unsigned int MeasurementVectorSize = 1;
    const unsigned int binsPerDimension = 30;

    HistogramType::MeasurementVectorType lowerBound(MeasurementVectorSize);
    lowerBound.Fill(0);

    HistogramType::MeasurementVectorType upperBound(MeasurementVectorSize);
    upperBound.Fill(255) ;

    HistogramType::SizeType size(MeasurementVectorSize);
    size.Fill(binsPerDimension);


    itk::ImageRegionConstIterator<FloatScalarImageType> imageIterator(indexSelectionFilter->GetOutput(), region);

    HistogramType::Pointer histogram = HistogramType::New();

    histogram->SetMeasurementVectorSize(MeasurementVectorSize);

    histogram->Initialize(size, lowerBound, upperBound );

    while(!imageIterator.IsAtEnd())
      {

      float pixel = imageIterator.Get();
      HistogramType::MeasurementVectorType mv(MeasurementVectorSize);
      mv[0] = pixel;

      histogram->IncreaseFrequencyOfMeasurement(mv, 1);

      ++imageIterator;
      }

    histograms.push_back(histogram);
    }

  return histograms;
}*/
/*
template <typename TBinValue>
HistogramType::Pointer Histogram<TBinValue>::ComputeNDHistogramOfRegionManual(const FloatVectorImageType* image, const itk::ImageRegion<2>& region, const unsigned int binsPerDimension)
{
  const unsigned int MeasurementVectorSize = image->GetNumberOfComponentsPerPixel();

  HistogramType::MeasurementVectorType lowerBound(MeasurementVectorSize);
  lowerBound.Fill(0);

  HistogramType::MeasurementVectorType upperBound(MeasurementVectorSize);
  upperBound.Fill(255) ;

  HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);

  itk::ImageRegionConstIterator<FloatVectorImageType> imageIterator(image, region);

  HistogramType::Pointer histogram = HistogramType::New();

  histogram->SetMeasurementVectorSize(MeasurementVectorSize);

  histogram->Initialize(size, lowerBound, upperBound );

  while(!imageIterator.IsAtEnd())
    {
    FloatVectorImageType::PixelType pixel = imageIterator.Get();
    HistogramType::MeasurementVectorType mv(MeasurementVectorSize);
    for(unsigned int i = 0; i < MeasurementVectorSize; ++i)
      {
      mv[i] = pixel[i];
      }

    histogram->IncreaseFrequencyOfMeasurement(mv, 1);

    ++imageIterator;
    }

  return histogram;
}*/
/*
template <typename TBinValue>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::Compute1DHistogramOfMultiChannelMaskedImage(const FloatVectorImageType* image, const itk::ImageRegion<2>& imageRegion, const Mask* mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins)
{
  // Compute the histogram for each channel separately
  std::vector<HistogramType::Pointer> channelHistograms = ComputeHistogramsOfMaskedRegion(image, imageRegion, mask, maskRegion, numberOfBins);

  unsigned int binsPerHistogram = channelHistograms[0]->GetSize(0);
  HistogramType histogram(binsPerHistogram*image->GetNumberOfComponentsPerPixel());

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
    {
    for(unsigned int bin = 0; bin < binsPerHistogram; ++bin)
      {
      histogram[channel*binsPerHistogram + bin] = channelHistograms[channel]->GetFrequency(bin);
      }
    }

  return histogram;
}*/
/*
template <typename TBinValue>
typename Histogram<TBinValue>::HistogramType Histogram<TBinValue>::Compute1DHistogramOfMultiChannelImage(const FloatVectorImageType* image, const itk::ImageRegion<2>& region, const unsigned int numberOfBins)
{
  // Compute the histogram for each channel separately
  std::vector<HistogramType::Pointer> channelHistograms = ComputeHistogramsOfRegionManual(image, region, numberOfBins);

  typename Histogram<TBinValue>::HistogramType histogram(channelHistograms[0]->GetSize(0));

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
    {
    for(unsigned int bin = 0; bin < channelHistograms[0]->GetSize(0); ++bin)
      {
      histogram.push_back(channelHistograms[channel]->GetFrequency(bin));
      }
    }

  return histogram;
}*/
/*
template <typename TBinValue>
HistogramType::Pointer Histogram<TBinValue>::ComputeNDHistogramOfMaskedRegionManual(const FloatVectorImageType* image, const Mask* mask, const itk::ImageRegion<2>& region, const unsigned int binsPerDimension)
{
  const unsigned int MeasurementVectorSize = image->GetNumberOfComponentsPerPixel();

  HistogramType::MeasurementVectorType lowerBound(MeasurementVectorSize);
  lowerBound.Fill(0);

  HistogramType::MeasurementVectorType upperBound(MeasurementVectorSize);
  upperBound.Fill(255) ;

  HistogramType::SizeType size(MeasurementVectorSize);
  size.Fill(binsPerDimension);

  itk::ImageRegionConstIterator<FloatVectorImageType> imageIterator(image, region);

  HistogramType::Pointer histogram = HistogramType::New();

  histogram->SetMeasurementVectorSize(MeasurementVectorSize);

  histogram->Initialize(size, lowerBound, upperBound );

  while(!imageIterator.IsAtEnd())
    {
    if(mask->IsHole(imageIterator.GetIndex()))
      {
      ++imageIterator;
      continue;
      }
    FloatVectorImageType::PixelType pixel = imageIterator.Get();
    HistogramType::MeasurementVectorType mv(MeasurementVectorSize);
    for(unsigned int i = 0; i < MeasurementVectorSize; ++i)
      {
      mv[i] = pixel[i];
      }

    histogram->IncreaseFrequencyOfMeasurement(mv, 1);

    ++imageIterator;
    }

  return histogram;
}*/
/*
template <typename TBinValue>
void Histogram<TBinValue>::OutputHistogram(const HistogramType::Pointer histogram)
{
  for(unsigned int i = 0; i < histogram->GetSize(0); ++i)
    {
    std::cout << histogram->GetFrequency(i) << " ";
    }
  std::cout << std::endl;
}*/
/*
template <typename TBinValue>
float Histogram<TBinValue>::NDHistogramDifference(const HistogramType::Pointer histogram1, const HistogramType::Pointer histogram2)
{
  unsigned int totalBins = 1;
  for(unsigned int i = 0; i < histogram1->GetSize().GetNumberOfElements(); ++i)
    {
    totalBins *= histogram1->GetSize(i);
    if(histogram1->GetSize(i) != histogram2->GetSize(i))
      {
      std::cerr << "Histograms must be the same size!" << std::endl;
      return 0;
      }
    }

  float totalDifference = 0;

  for(unsigned int i = 0; i < totalBins; ++i)
    {
    // The casts to float are necessary other wise the integer division always ends up = 0 !
    float normalized1 = static_cast<float>(histogram1->GetFrequency(i))/static_cast<float>(histogram1->GetTotalFrequency());
    float normalized2 = static_cast<float>(histogram2->GetFrequency(i))/static_cast<float>(histogram2->GetTotalFrequency());
    float difference = fabs(normalized1 - normalized2);
    //std::cout << "difference: " << difference << std::endl;
    totalDifference += difference;
    }

  //std::cout << "totalDifference: " << totalDifference << std::endl;
  return totalDifference;
}*/
/*
template <typename TBinValue>
float Histogram<TBinValue>::HistogramDifference(const HistogramType::Pointer histogram1, const HistogramType::Pointer histogram2)
{
  if(histogram1->GetSize(0) != histogram2->GetSize(0))
    {
    std::cerr << "Histograms must be the same size!" << std::endl;
    return 0;
    }

  float totalDifference = 0;
  for(unsigned int i = 0; i < histogram1->GetSize(0); ++i)
    {
//     std::cout << "h1 " << i << " " << histogram1->GetFrequency(i) << std::endl;
//     std::cout << "h1 total " << histogram1->GetTotalFrequency() << std::endl;
//
//     std::cout << "h2 " << i << " " << histogram2->GetFrequency(i) << std::endl;
//     std::cout << "h2 total " << histogram2->GetTotalFrequency() << std::endl;

    // The casts to float are necessary other wise the integer division always ends up = 0 !
    float normalized1 = static_cast<float>(histogram1->GetFrequency(i))/static_cast<float>(histogram1->GetTotalFrequency());
    float normalized2 = static_cast<float>(histogram2->GetFrequency(i))/static_cast<float>(histogram2->GetTotalFrequency());
    //float difference = fabs(normalized1 - normalized2);
    float difference = (normalized1 - normalized2)*(normalized1 - normalized2);
    //std::cout << "difference: " << difference << std::endl;
    totalDifference += difference;
    }

  //std::cout << "totalDifference: " << totalDifference << std::endl;
  return totalDifference;
}*/
/*
template <typename TBinValue>
float Histogram<TBinValue>::HistogramIntersection(const HistogramType::Pointer histogram1, const HistogramType::Pointer histogram2)
{
  if(histogram1->GetSize(0) != histogram2->GetSize(0))
    {
    std::cerr << "Histograms must be the same size!" << std::endl;
    return 0;
    }

  float totalIntersection = 0.0f;
  for(unsigned int i = 0; i < histogram1->GetSize(0); ++i)
    {
    // The casts to float are necessary other wise the integer division always ends up = 0 !
    float frequency1 = static_cast<float>(histogram1->GetFrequency(i));
    float frequency2 = static_cast<float>(histogram2->GetFrequency(i));
    //float difference = fabs(normalized1 - normalized2);
    float interseciton = std::min(frequency1, frequency2);
    //std::cout << "difference: " << difference << std::endl;
    totalIntersection += interseciton;
    }

  float normalizedIntersection = totalIntersection / static_cast<float>(histogram2->GetTotalFrequency());

  return normalizedIntersection;
}*/



#endif
