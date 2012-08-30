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

#ifndef MaskedHistogram_HPP
#define MaskedHistogram_HPP

#include "MaskedHistogram.h"

template <typename TBinValue>
template <typename TComponent>
typename MaskedHistogramGenerator<TBinValue>::HistogramType MaskedHistogramGenerator<TBinValue>::ComputeMaskedImage1DHistogram(
                const itk::VectorImage<TComponent, 2>* image,
                const itk::ImageRegion<2>& imageRegion,
                const Mask* const mask,
                const itk::ImageRegion<2>& maskRegion,
                const unsigned int numberOfBinsPerDimension,
                const TComponent& rangeMin,
                const TComponent& rangeMax, const bool allowOutside, const unsigned char maskValue)
{
  // For VectorImage, we must use VectorImageToImageAdaptor

  typedef itk::VectorImage<TComponent, 2> ImageType;
  typedef itk::Image<TComponent, 2> ScalarImageType;

  std::vector<itk::Index<2> > maskIndices = ITKHelpers::GetPixelsWithValueInRegion(mask, maskRegion, maskValue);

  // Compute the corresponding locations in the imageRegion
  std::vector<itk::Index<2> > imageIndices(maskIndices.size());
  itk::Offset<2> maskRegionToImageRegionOffset = imageRegion.GetIndex() - maskRegion.GetIndex(); // The offset from the maskRegion to the imageRegion
  for(size_t i = 0; i < maskIndices.size(); ++i)
  {
    imageIndices[i] = maskIndices[i] + maskRegionToImageRegionOffset;
  }

  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel

    typedef itk::VectorImageToImageAdaptor<TComponent, 2> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SetExtractComponentIndex(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    // Get this channels masked scalar values
    std::vector<typename ScalarImageType::PixelType> validPixels =
        ITKHelpers::GetPixelValues(adaptor.GetPointer(), imageIndices);

    // Compute the histogram of the scalar values
    HistogramType histogram = HistogramGeneratorType::ScalarHistogram(validPixels, numberOfBinsPerDimension, rangeMin, rangeMax, allowOutside);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue>
template <typename TComponent, unsigned int Dimension>
typename MaskedHistogramGenerator<TBinValue>::HistogramType MaskedHistogramGenerator<TBinValue>::ComputeMaskedImage1DHistogram
    (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
     const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBinsPerDimension,
     const TComponent& rangeMin, const TComponent& rangeMax, const bool allowOutside, const unsigned char maskValue)
{
  // For Image<CovariantVector>, we must use NthElementImageAdaptor

  typedef itk::Image<itk::CovariantVector<TComponent, Dimension>, 2> ImageType;
  typedef itk::Image<TComponent, 2> ScalarImageType;

  std::vector<itk::Index<2> > maskIndices = ITKHelpers::GetPixelsWithValueInRegion(mask, maskRegion, maskValue);

  // Compute the corresponding locations in the imageRegion
  std::vector<itk::Index<2> > imageIndices(maskIndices.size());
  itk::Offset<2> maskRegionToImageRegionOffset = imageRegion.GetIndex() - maskRegion.GetIndex(); // The offset from the maskRegion to the imageRegion
  for(size_t i = 0; i < maskIndices.size(); ++i)
  {
    imageIndices[i] = maskIndices[i] + maskRegionToImageRegionOffset;
  }

  HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < Dimension; ++channel)
  {
    // Extract the channel

    typedef itk::NthElementImageAdaptor<ImageType, typename ScalarImageType::PixelType> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SelectNthElement(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    // Get this channels masked scalar values
    std::vector<typename ScalarImageType::PixelType> validPixels =
        ITKHelpers::GetPixelValues(adaptor.GetPointer(), imageIndices);

    HistogramType histogram;
    // Compute the histogram of the scalar values
    if(allowOutside)
    {
      histogram = HistogramGeneratorType::ScalarHistogramAllowOutside(validPixels, numberOfBinsPerDimension, rangeMin, rangeMax);
    }
    else
    {
      histogram = HistogramGeneratorType::ScalarHistogram(validPixels, numberOfBinsPerDimension, rangeMin, rangeMax);
    }

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TBinValue>
template <typename TComponent, unsigned int Dimension>
typename MaskedHistogramGenerator<TBinValue>::HistogramType MaskedHistogramGenerator<TBinValue>::ComputeQuadrantMaskedImage1DHistogram
    (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
     const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBinsPerDimension,
     const TComponent& rangeMin, const TComponent& rangeMax, const bool allowOutside, const unsigned char maskValue)
{
  HistogramType fullHistogram;

  // Ignore the quadrant if it doesn't contain at least this ratio of valid pixels (i.e. too much of the quadrant is the hole)
  float minValidPixelRatio = 0.2f; // If less than 20% of the pixels in a quadrant are valid, do not consider the quadrant at all

  for(unsigned int quadrant = 0; quadrant < 4; ++quadrant)
  {
    itk::ImageRegion<2> maskRegionQuadrant = ITKHelpers::GetQuadrant(maskRegion, quadrant);
    itk::ImageRegion<2> imageRegionQuadrant = ITKHelpers::GetQuadrant(imageRegion, quadrant);

    // Compute the ratio of valid pixels to total pixels in the quadrant
    float validPixelRatio = static_cast<float>(mask->CountValidPixels(maskRegionQuadrant)) / static_cast<float>(maskRegionQuadrant.GetNumberOfPixels());

    // Only calculate and append the histogram if there are enough valid pixels
    if(validPixelRatio > minValidPixelRatio)
    {
      HistogramType quadrantHistogram = ComputeMaskedImage1DHistogram(image, imageRegionQuadrant, mask, maskRegionQuadrant, numberOfBinsPerDimension, rangeMin, rangeMax, allowOutside, maskValue);
      fullHistogram.insert(fullHistogram.end(), quadrantHistogram.begin(), quadrantHistogram.end());
    }
  }

//  std::cout << "Histogram has " << fullHistogram.size() << " bins." << std::endl;
  return fullHistogram;
}

#endif
