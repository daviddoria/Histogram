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

#ifndef QuadrantHistogram_H
#define QuadrantHistogram_H

#include "Histogram.h"

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

  unsigned int NumberOfBinsPerDimension;

  /** This flag determines if an error is produced if a value is found to be outside the requested range.
    * If it is set to "false", the error is produced. If it is set to "true", values outside the range are
    * added to the closesest end bin (highest bin or lowest bin, depending on if the value is above or below
    * the requested range, respectively). */
  bool AllowOutside;


  /** A boolean for each quadrant indicating if it was computed (if the image patch had enough valid pixels to be used). */
  bool Valid[4];

  QuadrantHistogramProperties() : NumberOfBinsPerDimension(30), AllowOutside(true)
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      this->Valid[i] = false;
    }
  }
};

/** This class stores 4 histograms, representing the histograms of the 4 quadrants of a region. */
template <typename THistogram, typename TQuadrantProperties>
struct QuadrantHistogram
{
  TQuadrantProperties Properties;

  THistogram Histograms[4];

  QuadrantHistogram()
  {

  }

  void NormalizeHistograms()
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      if(this->Properties.Valid[i])
      {
        Helpers::NormalizeVectorInPlace(this->Histograms[i]);
      }
    }
  }
};

#endif
