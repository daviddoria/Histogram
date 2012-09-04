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

/** This class stores properties about quadrant histograms. */
template <TPixelType>
struct QuadrantHistogramProperties
{
  typedef std::vector<TPixelType> RangeContainerType;

  /** A collection of 4 min ranges (one for each quadrant). */
  std::vector<RangeContainerType> QuadrantMinRanges;

  /** A collection of 4 max ranges (one for each quadrant). */
  std::vector<RangeContainerType> QuadrantMaxRanges;

  /** A boolean for each quadrant indicating if it has enough valid pixels to be used. */
  std::vector<bool> Used;
};

/** This collection of functions computes histograms for regions of images where only some of the pixels
  * are considered, based on their corresponding value in a Mask. */
template <typename TBinValue>
class MaskedHistogramGenerator
{
public:

  /** The type of histogram generator to use to generate non-masked histograms internally. */
  typedef HistogramGenerator<TBinValue> HistogramGeneratorType;

  /** The type of histogram that will be returned by functions in this class. */
  typedef typename HistogramGeneratorType::HistogramType HistogramType;

  /** This function computes the histogram of valid/hole (specified by 'maskValue') pixels (according to 'mask') in an image.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::VectorImage. */
  template <typename TComponent>
  static HistogramType ComputeMaskedImage1DHistogram
      (const itk::VectorImage<TComponent, 2>* const image,
       const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask,
       const itk::ImageRegion<2>& maskRegion,
       const unsigned int numberOfBins,
       const std::vector<TComponent>& rangeMins,
       const std::vector<TComponent>& rangeMaxs,
       const bool allowOutside,
       const unsigned char maskValue);

  /** This function computes the histogram of valid/hole (specified by 'maskValue') pixels (according to 'mask') in an image.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>. */
  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeMaskedImage1DHistogram
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
       const std::vector<TComponent>& rangeMins, const std::vector<TComponent>& rangeMaxs, const bool allowOutside, const unsigned char maskValue);

  /** This function computes and concatenates the histogram of valid/hole (specified by 'maskValue') pixels
    * (according to 'mask') in an image in each of the 4 quadrants, using the same specified range for each channel.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>.
    * For example, if a 3 channel image is presented and 10 bins per channel are requested, the output histogram will be
    * 4 (quadrants) * 3 (channels) * 10 (bins per channel) = 120 bins. */
  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeQuadrantMaskedImage1DHistogram
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
       const TComponent& rangeMin, const TComponent& rangeMax, const bool allowOutside, const unsigned char maskValue);

  /** This function computes and concatenates the histogram of valid/hole (specified by 'maskValue') pixels
    * (according to 'mask') in an image in each of the 4 quadrants using a different specified range for each channel.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>.
    * For example, if a 3 channel image is presented and 10 bins per channel are requested, the output histogram will be
    * 4 (quadrants) * 3 (channels) * 10 (bins per channel) = 120 bins.
    * TComponent and Dimension specify the template parameters of the itk::CovariantVector pixel type.
    * TRangeContainer is a vector container that must have an operator[]()
    */
  template <typename TComponent, unsigned int Dimension, typename TRangeContainer>
  static HistogramType ComputeQuadrantMaskedImage1DHistogram
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
       const TRangeContainer& rangeMins, const TRangeContainer& rangeMaxs,
       const bool allowOutside, const unsigned char maskValue);

  /** This function computes and concatenates the histogram of valid/hole (specified by 'maskValue') pixels
    * (according to 'mask') in an image in each of the 4 quadrants using an adaptive range for each channel of each quadrant.
    * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
    * a target mask to a source patch. This function is for itk::Image<CovariantVector>.
    * For example, if a 3 channel image is presented and 10 bins per channel are requested, the output histogram will be
    * 4 (quadrants) * 3 (channels) * 10 (bins per channel) = 120 bins. */
  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeQuadrantMaskedImage1DHistogramAdaptive
      (const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* const image, const itk::ImageRegion<2>& imageRegion,
       const Mask* const mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
       const bool allowOutside, const unsigned char maskValue);

}; // end class MaskedHistogramGenerator

#include "MaskedHistogramGenerator.hpp"

#endif
