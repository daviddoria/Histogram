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

#ifndef MaskedHistogram_H
#define MaskedHistogram_H

#include "MaskedHistogram.h"

namespace MaskedHistogram
{

template <typename TComponent>
Histogram<int>::HistogramType ComputeMaskedImageHistogram1D(
                const itk::VectorImage<TComponent, 2>* image,
                const itk::ImageRegion<2>& region,
                const Mask* const mask, const itk::ImageRegion<2>& maskRegion,
                const unsigned int numberOfBinsPerDimension,
                const TComponent& rangeMin,
                const TComponent& rangeMax)
{
  // For VectorImage, we must use VectorImageToImageAdaptor
  std::vector<itk::Index<2> > validIndices = mask->GetValidPixelsInRegion(maskRegion);

  Histogram<int>::HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel
    typedef itk::Image<typename TypeTraits<typename ImageType::PixelType>::ComponentType, 2> ScalarImageType;

    typedef itk::VectorImageToImageAdaptor<typename TypeTraits<typename ImageType::PixelType>::ComponentType, 2> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SetExtractComponentIndex(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    // Get this channels masked scalar values
    std::vector<typename ScalarImageType::PixelType> validPixels =
        ITKHelpers::GetPixelValues(adaptor.GetPointer(), validIndices);

    // Compute the histogram of the scalar values
    HistogramType histogram = ScalarHistogram(validPixels, numberOfBinsPerDimensions, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TImage, typename TMask>
Histogram<int> ComputeMaskedImage1DHistogram(const TImage* const image, const itk::ImageRegion<2>& imageRegion,
                                               const Mask* const mask, const itk::ImageRegion<2>& maskRegion,
                                               const unsigned int numberOfBins, const typename TImage::PixelType& rangeMin, const typename TImage::PixelType& rangeMax)
{
  throw std::runtime_error("ComputeMaskedImage1DHistogram is not yet implemented for generic image types.!")
}

} // end namespace

#endif
