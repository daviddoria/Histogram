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

#ifndef MaskedHistogramGenerator_H
#define MaskedHistogramGenerator_H

#include "HistogramGenerator.h"

// STL
#include <vector>

/** This class stores properties about quadrant histograms.
  * \tparam TRangeContainer - A type that holds a scalar value for each dimension of the image.
  * This is typically TImage::PixelType or std::vector<typename TypeTraits<typename TImage::PixelType>::ComponentType>.
  */
template <typename TRangeContainer>
struct QuadrantHistogramProperties
{
  /** A collection of min ranges (one for each quadrant). */
  TRangeContainer QuadrantMinRanges[4];

  /** A collection of max ranges (one for each quadrant). */
  TRangeContainer QuadrantMaxRanges[4];

  /** A boolean for each quadrant indicating if it has enough valid pixels to be used. */
  bool Used[4];

  unsigned int NumberOfBinsPerDimension;

  /** This flag determines if an error is produced if a value is found to be outside the requested range.
    * If it is set to "false", the error is produced. If it is set to "true", values outside the range are
    * added to the closesest end bin (highest bin or lowest bin, depending on if the value is above or below
    * the requested range, respectively). */
  bool AllowOutside;

  QuadrantHistogramProperties() : NumberOfBinsPerDimension(30), AllowOutside(true)
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      this->Used[i] = false;
    }
  }
};

/** This collection of functions computes histograms for regions of images where only some of the pixels
  * are considered, based on their corresponding value in a Mask. */
template <typename TBinValue>
class MaskedHistogramGenerator
{
public:

  /** The type of histogram generator to use to generate non-masked histograms internally. */
  typedef HistogramGenerator<TBinValue> HistogramGeneratorType;

  /** The type of histogram that will be returned by single histogram functions in this class. */
  typedef typename HistogramGeneratorType::HistogramType HistogramType;

  /** The type of histogram that will be returned by quadrant histogram functions in this class. */
  typedef typename HistogramGeneratorType::QuadrantHistogramType QuadrantHistogramType;

  /** This function computes the histogram of valid/hole (specified by 'maskValue') pixels (according to 'mask') in an image.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::VectorImage. */
  template <typename TComponent, typename TRangeContainer>
  static HistogramType ComputeMaskedImage1DHistogram
      (const itk::VectorImage<TComponent, 2>* const image,
       const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask,
       const itk::ImageRegion<2>& maskRegion,
       const unsigned int numberOfBins,
       const TRangeContainer& rangeMins,
       const TRangeContainer& rangeMaxs,
       const bool allowOutside,
       const unsigned char maskValue);

  /** This function computes the histogram of valid/hole (specified by 'maskValue') pixels (according to 'mask') in an image.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>. */
  template <typename TComponent, unsigned int Dimension, typename TRangeContainer>
  static HistogramType ComputeMaskedImage1DHistogram
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
       const TRangeContainer& rangeMins, const TRangeContainer& rangeMaxs, const bool allowOutside, const unsigned char maskValue);

  /** This function computes and concatenates the histogram of valid/hole (specified by 'maskValue') pixels
    * (according to 'mask') in an image in each of the 4 quadrants using an adaptive range for each channel of each quadrant.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>.
    * For example, if a 3 channel image is presented and 10 bins per channel are requested, the output histogram will be
    * 4 (quadrants) * 3 (channels) * 10 (bins per channel) = 120 bins.
    *
    * We return the details about the histogram in returnQuadrantHistogramProperties so that they can be used to compute
    * future corresponding histograms.
    *
    * 'useProvidedRanges' determines if the ranges (and 'Used' flags) that are in the quadrantHistogramProperties are used,
    * or if the ranges are computed adaptively. The 'NumberOfComponentsPerDimension' and 'AllowOutside' are always used.
    */
  template <typename TComponent, unsigned int Dimension>
  static QuadrantHistogramType ComputeQuadrantMaskedImage1DHistogramAdaptive
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const QuadrantHistogramProperties<itk::CovariantVector<TComponent, Dimension> >& quadrantHistogramProperties,
       const bool useProvidedRanges,
       const unsigned char maskValue, QuadrantHistogramProperties<itk::CovariantVector<TComponent, Dimension> >& returnQuadrantHistogramProperties);

}; // end class MaskedHistogramGenerator

#include "MaskedHistogramGenerator.hpp"

#endif
