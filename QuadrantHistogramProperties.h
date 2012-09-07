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

#ifndef QuadrantHistogramProperties_H
#define QuadrantHistogramProperties_H

/** This class stores properties about quadrant histograms.
  * \tparam TRangeContainer - A type that holds a scalar value for each dimension of the image.
  * This is typically TImage::PixelType or std::vector<typename TypeTraits<typename TImage::PixelType>::ComponentType>.
  */
template <typename TRangeContainer = std::vector<float> >
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

  /** This flag lets the user know if there has been any previous activity in this instance. This lets us compute
    * which quadrants have enough valid pixels the first time this set of properties is used, etc. That is, if we specify
    * the ranges for all of the quadrants, we could not just set all of the quadrants to Valid, because they might not have
    * enough pixels! */
  bool Initialized;

  /** A boolean for each quadrant indicating if it was computed (if the image patch had enough valid pixels to be used). */
  bool Valid[4];

  QuadrantHistogramProperties() : NumberOfBinsPerDimension(30), AllowOutside(true), Initialized(false)
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      this->Valid[i] = false;
    }
  }

  void SetAllMinRanges(const TRangeContainer& minRange)
  {
    for(unsigned int quadrant = 0; quadrant < 4; ++quadrant)
    {
      this->QuadrantMinRanges[quadrant] = minRange;
    }
  }

  void SetAllMaxRanges(const TRangeContainer& maxRange)
  {
    for(unsigned int quadrant = 0; quadrant < 4; ++quadrant)
    {
      this->QuadrantMaxRanges[quadrant] = maxRange;
    }
  }
};

#endif
