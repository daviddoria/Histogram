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

#ifndef HistogramGenerator_H
#define HistogramGenerator_H

// STL
#include <string>
#include <vector>

// ITK
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkVectorImage.h"

// Submodules
#include "Helpers/TypeTraits.h"

#include "Histogram.h"
#include "QuadrantHistogram.h"

/** A collection of functions to compute histograms.
  * These are static functions in a class so that the bin type of the returned histograms can be specified more easily (as the class template paramter).
  *
  * This class is named HistogramGenerator to avoid naming confusion between HistogramType and HistogramGeneratorType.
  * That is, in user code we will have:
  * typedef HistogramGenerator<int> HistogramGeneratorType;
  * typedef HistogramGeneratorType::HistogramType HistogramType;
  *
  * which makes sense, versus:
  * typedef Histogram<int> HistogramType;
  * typedef HistogramGeneratorType::HistogramType [what would this be called?];
  *
  *  \tparam TQuadrantProperties the type of QuadrantHistogramProperties to use. A default type is provided so that
  *  in cases where we are not using quadrants, no type must be specified.
  */
template <typename TBinValue, typename TQuadrantProperties = std::vector<float> >
class HistogramGenerator
{
public:
  typedef TBinValue BinValueType;

//  typedef std::vector<TBinValue> HistogramType;
  typedef Histogram<TBinValue> HistogramType;

  typedef QuadrantHistogram<HistogramType, TQuadrantProperties> QuadrantHistogramType;

  /** Compute the histogram of a scalar itk::Image. This function is called by the ComputeImageHistogram1D overload
    * that does not match itk::VectorImage or itk::Image<CovariantVector>. */
  template <typename TImage>
  static HistogramType ComputeScalarImageHistogram(
                  const TImage* image, const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
                  const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax, const bool allowOutside = false);

  /** Compute the histogram of a scalar itk::Image. */
  template <typename TScalarImage>
  static HistogramType ComputeImageHistogram1D(
                  const TScalarImage* image, const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMin,
                  const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMax, const bool allowOutside = false);

  /** Compute the histograms of each channel of an itk::Image<CovariantVector>, and concatentate them together to form
    * of a 1D histogram. It is necessary to have this function separate from the one that handles itk::VectorImage because
    * NthElementImageAdaptor must be used for Image<CovariantVector>, while VectorImageToImageAdapter must be used for VectorImage. */
  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeImageHistogram1D(
                  const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* image,
                  const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const TComponent& rangeMin,
                  const TComponent& rangeMax, const bool allowOutside = false);

  /** Compute the histograms of each channel of an itk::VectorImage, and concatentate them together to form
    * of a 1D histogram.  It is necessary to have this function separate from the one that handles itk::Image<CovariantVector> because
    * NthElementImageAdaptor must be used for Image<CovariantVector>, while VectorImageToImageAdapter must be used for VectorImage. */
  template <typename TComponent>
  static HistogramType ComputeImageHistogram1D(
                  const itk::VectorImage<TComponent, 2>* image,
                  const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const TComponent& rangeMin,
                  const TComponent& rangeMax, const bool allowOutside);

  /** Compute the histogram of a collection of values. The values are expected to fall inside of [rangeMin, rangeMax].*/
  template <typename TValue>
  static HistogramType ScalarHistogram(const std::vector<TValue>& values, const unsigned int numberOfBins,
                                const TValue& rangeMin, const TValue& rangeMax);

  /** Compute the histogram of a collection of values. Values that fall outside of [rangeMin, rangeMax] are placed in the end bin
    * closest to them. That is, if the range is [0,10] and the value is 11, it will be placed in the bin in which the top of the bin range
    * is 10. Conversely, if the range is [0, 10] and the value is -2, it will be placed in the bin in which the bottom of the bin range is 0.*/
  template <typename TValue>
  static HistogramType ScalarHistogramAllowOutside(const std::vector<TValue>& values, const unsigned int numberOfBins,
                                const TValue& rangeMin, const TValue& rangeMax);

  /** Compute the joint histogram of the channels of an itk::Image<CovariantVector>. Though the histogram is ND (where N is the number of channels),
    * we output a 1D histogram (the ND histgram linearized in a "raster scan" order). If a 'mask' is provided, the pixels with corresponding pixels
    * in 'mask' equal to 'maskValue' are the only ones used. */
  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeImageJointChannelHistogram(
      const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* image,
      const itk::ImageRegion<2>& region,
      const unsigned int numberOfBinsPerDimension,
      const itk::CovariantVector<TComponent, Dimension>& rangeMin, const itk::CovariantVector<TComponent, Dimension>& rangeMax,
      const Mask* const mask = 0, const Mask::PixelType& maskValue = 255,
      const itk::ImageRegion<2>& maskRegion = itk::ImageRegion<2>(),
      const bool allowOutside = false);

}; // end class HistogramGenerator

#include "HistogramGenerator.hpp"

#endif
