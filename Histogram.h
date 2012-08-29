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

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

// STL
#include <string>
#include <vector>

// ITK
#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkVectorImage.h"

// Submodules
#include "Helpers/TypeTraits.h"

template <typename TBinValue>
class Histogram
{
public:
  typedef std::vector<TBinValue> HistogramType;

  /** Compute the histograms of each channel of an image, and concatentate them together to form
    * of a 1D histogram. */

  template <typename TImage>
  static HistogramType ComputeScalarImageHistogram(
                  const TImage* image, const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMin,
                  const typename TypeTraits<typename TImage::PixelType>::ComponentType& rangeMax);
  
  template <typename TScalarImage>
  static HistogramType ComputeImageHistogram1D(
                  const TScalarImage* image, const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMin,
                  const typename TypeTraits<typename TScalarImage::PixelType>::ComponentType& rangeMax);

  template <typename TComponent, unsigned int Dimension>
  static HistogramType ComputeImageHistogram1D(
                  const itk::Image<itk::CovariantVector<TComponent, Dimension>, 2>* image,
                  const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const TComponent& rangeMin,
                  const TComponent& rangeMax);

  template <typename TComponent>
  static HistogramType ComputeImageHistogram1D(
                  const itk::VectorImage<TComponent, 2>* image,
                  const itk::ImageRegion<2>& region,
                  const unsigned int numberOfBinsPerDimension,
                  const TComponent& rangeMin,
                  const TComponent& rangeMax);

  /** Compute the histogram of a collection of values. */
  template <typename TValue>
  static HistogramType ScalarHistogram(const std::vector<TValue>& values, const unsigned int numberOfBins,
                                const TValue& rangeMin, const TValue& rangeMax);

  /** Compute the "Histogram Intersection" score between two histograms. This is the sum of the minimum of each corresponding bin.
    * There is also a normalization that is performed.
    * (M. J. Swain and D. H. Ballard. "Color Indexing." IJCV 7(1):11-32 November 1991) */
  static float HistogramIntersection(const HistogramType& histogram1, const HistogramType& histogram2);

  /** Compute a bin-to-bin difference between the histograms. */
  static float HistogramDifference(const HistogramType& histogram1, const HistogramType& histogram2);

  /** Compute a weighted bin-to-bin difference between the histograms. Bins with large counts in 'idealHistogram'
    * are considered to have low weight, as the difference at this bin is not necessarily bad (just because a certain color doesn't
    * appear in a query histogram, that doesn't mean it is a bad match.) However, differences at bins with small counts in 'idealHistogram'
    * are important (contribute a lot to the error) because this means that new colors have been introduced, which is bad for "color coherence." */
  static float WeightedHistogramDifference(const HistogramType& idealHistogram, const HistogramType& queryHistogram);

  /** Determine how many bins in 'queryHistogram' have non-zero counts while the corresponding bin in 'idealHistogram' has a zero count. The idea is to get a sense of
    * how many colors are present in the queryHistogram that are not present in the idealHistogram. A low score means that all of the colors in the queryHistogram are
    * present in the idealHistogram, while a high score means that different colors are present. */
  static float HistogramCoherence(const HistogramType& idealHistogram, const HistogramType& queryHistogram);

  /** Write the histogram bin counts to a file. */
  static void WriteHistogram(const HistogramType& histogram, const std::string& filename);

  /** Write the histogram bin counts to the screen. */
  static void OutputHistogram(const HistogramType& histogram);

};

#include "Histogram.hpp"

#endif
