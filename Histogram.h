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

#ifndef Histogram_H
#define Histogram_H

// Submodules
#include <Helpers/Helpers.h>

// STL
#include <vector>

/** This class is a simple wrapper around std::vector to give it functions
  * named appropriately for histogram storage. */
template <typename TBinValue>
class Histogram : public std::vector<TBinValue>
{
public:

  typedef TBinValue BinValueType;

  typedef std::vector<TBinValue> VectorType;

  Histogram() : VectorType(){}

  Histogram(const unsigned int numberOfBins) : VectorType(numberOfBins){}

  Histogram(const unsigned int numberOfBins, const TBinValue binValue) : VectorType(numberOfBins, binValue){}

  void Normalize()
  {
    Helpers::NormalizeVectorInPlace(*this);
  }

  void Print()
  {
    Helpers::Output(*this);
  }

  void PrintPadded()
  {
    for(size_t i = 0; i < this->size(); ++i)
    {
      std::cout << Helpers::ZeroPad(this->at(i), 3);
    }
  }

};

#endif
