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

#include "Histogram.h"

namespace MaskedHistogram
{
/** This function computes the histogram of valid pixels (according to 'mask') in an image.
  * The 'maskRegion' is not necessarily the same as the 'imageRegion', as we may want to apply
  * a target mask to a source patch. */
template <typename TImage, typename TMask>
Histogram<int>::HistogramType ComputeMaskedImage1DHistogram
    (const TImage* image, const itk::ImageRegion<2>& imageRegion,
     Mask* mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins,
     const typename TImage::PixelType& rangeMin, const typename TImage::PixelType& rangeMax);

}

#endif
