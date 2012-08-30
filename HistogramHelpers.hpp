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

#ifndef HistogramHelpers_H
#define HistogramHelpers_H

// STL
#include <string>
#include <vector>

namespace HistogramHelpers
{
  /** Write the histogram bin counts to a file. */
  template <typename THistogram>
  static void WriteHistogram(const THistogram& histogram, const std::string& filename)
  {
    std::ofstream fout(filename.c_str());
    for(unsigned int i = 0; i < histogram.size(); ++i)
    {
      fout << histogram[i] << " ";
    }

    fout.close();
  }

  /** Write the histogram bin counts to the screen. */
  template <typename THistogram>
  static void OutputHistogram(const THistogram& histogram)
  {
    for(unsigned int i = 0; i < histogram.size(); ++i)
    {
      std::cout << histogram[i] << " ";
    }
  }
}

#endif
