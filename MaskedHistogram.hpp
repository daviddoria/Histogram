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

namespace MaskedHistogram
{

template <typename TComponent>
Histogram<int>::HistogramType ComputeMaskedImageHistogram1D(
                const itk::VectorImage<TComponent, 2>* image,
                const itk::ImageRegion<2>& imageRegion,
                const Mask* const mask, const itk::ImageRegion<2>& maskRegion,
                const unsigned int numberOfBinsPerDimension,
                const TComponent& rangeMin,
                const TComponent& rangeMax)
{
  // For VectorImage, we must use VectorImageToImageAdaptor

  typedef itk::VectorImage<TComponent, 2> ImageType;
  typedef itk::Image<TComponent, 2> ScalarImageType;

  std::vector<itk::Index<2> > validMaskIndices = mask->GetValidPixelsInRegion(maskRegion);

  // Compute the corresponding locations in the imageRegion
  std::vector<itk::Index<2> > validImageIndices(validMaskIndices.size());
  itk::Offset<2> maskRegionToImageRegionOffset = imageRegion.GetIndex() - maskRegion.GetIndex(); // The offset from the maskRegion to the imageRegion
  for(size_t i = 0; i < validMaskIndices.size(); ++i)
  {
    validImageIndices[i] = validMaskIndices[i] + maskRegionToImageRegionOffset;
  }

  Histogram<int>::HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
  {
    // Extract the channel

    typedef itk::VectorImageToImageAdaptor<TComponent, 2> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SetExtractComponentIndex(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    // Get this channels masked scalar values
    std::vector<typename ScalarImageType::PixelType> validPixels =
        ITKHelpers::GetPixelValues(adaptor.GetPointer(), validImageIndices);

    // Compute the histogram of the scalar values
    Histogram<int>::HistogramType histogram = Histogram<int>::ScalarHistogram(validPixels, numberOfBinsPerDimension, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

template <typename TComponent, unsigned int Dimension>
Histogram<int>::HistogramType ComputeMaskedImage1DHistogram
    (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
     const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBinsPerDimension,
     const TComponent& rangeMin, const TComponent& rangeMax)
{
  // For Image<CovariantVector>, we must use NthElementImageAdaptor

  typedef itk::Image<itk::CovariantVector<TComponent, Dimension>, 2> ImageType;
  typedef itk::Image<TComponent, 2> ScalarImageType;

  std::vector<itk::Index<2> > validMaskIndices = mask->GetValidPixelsInRegion(maskRegion);

  // Compute the corresponding locations in the imageRegion
  std::vector<itk::Index<2> > validImageIndices(validMaskIndices.size());
  itk::Offset<2> maskRegionToImageRegionOffset = imageRegion.GetIndex() - maskRegion.GetIndex(); // The offset from the maskRegion to the imageRegion
  for(size_t i = 0; i < validMaskIndices.size(); ++i)
  {
    validImageIndices[i] = validMaskIndices[i] + maskRegionToImageRegionOffset;
  }

  Histogram<int>::HistogramType concatenatedHistograms;

  for(unsigned int channel = 0; channel < Dimension; ++channel)
  {
    // Extract the channel

    typedef itk::NthElementImageAdaptor<ImageType, typename ScalarImageType::PixelType> ImageAdaptorType;
    typename ImageAdaptorType::Pointer adaptor = ImageAdaptorType::New();
    adaptor->SelectNthElement(channel);
    adaptor->SetImage(const_cast<ImageType*>(image));

    // Get this channels masked scalar values
    std::vector<typename ScalarImageType::PixelType> validPixels =
        ITKHelpers::GetPixelValues(adaptor.GetPointer(), validImageIndices);

    // Compute the histogram of the scalar values
    Histogram<int>::HistogramType histogram = Histogram<int>::ScalarHistogram(validPixels, numberOfBinsPerDimension, rangeMin, rangeMax);

    concatenatedHistograms.insert(concatenatedHistograms.end(), histogram.begin(), histogram.end());
  }

  return concatenatedHistograms;
}

} // end namespace

#endif
