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

#include "QuadrantHistogramProperties.h"

/** This class stores 4 histograms, representing the histograms of the 4 quadrants of a region. */
template <typename THistogram, typename TQuadrantProperties>
struct QuadrantHistogram
{
  typedef THistogram HistogramType;
  typedef TQuadrantProperties QuadrantPropertiesType;

  /** Store the properties that were used to compute the histograms. */
  TQuadrantProperties Properties;

  /** Store the 4 histograms. */
  THistogram Histograms[4];

  /** Constructor. */
  QuadrantHistogram()
  {

  }

  /** Normalize all of the valid histograms. */
  void NormalizeHistograms()
  {
    for(unsigned int i = 0; i < 4; ++i)
    {
      if(this->Properties.Valid[i])
      {
        this->Histograms[i].Normalize();
      }
    }
  }

  /** Print all of the valid histograms to stdout. */
  void PrintHistograms()
  {
    std::cout << "Valid quadrants: ";
    for(unsigned int i = 0; i < 4; ++i)
    {
      if(this->Properties.Valid[i])
      {
        std::cout << i << " ";
      }
    }

    std::cout << std::endl;

    for(unsigned int i = 0; i < 4; ++i)
    {
      if(this->Properties.Valid[i])
      {
        std::cout << "Quadrant " << i << ": ";
        this->Histograms[i].Print();
      }
    }
  }

};

#endif
