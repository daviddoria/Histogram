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

#ifndef HistogramDifferences_H
#define HistogramDifferences_H

// STL
#include <string>
#include <vector>

namespace HistogramDifferences
{

//  /** Compute the "Histogram Intersection" score between two histograms. This is the sum of the minimum of each corresponding bin.
//    * There is also a normalization that is performed.
//    * (M. J. Swain and D. H. Ballard. "Color Indexing." IJCV 7(1):11-32 November 1991) */
//  static float HistogramIntersection(const HistogramType& histogram1, const HistogramType& histogram2);

//  /** Compute a bin-to-bin difference between the histograms. */
//  static float HistogramDifference(const HistogramType& histogram1, const HistogramType& histogram2);

//  /** Compute a weighted bin-to-bin difference between the histograms. Bins with large counts in 'idealHistogram'
//    * are considered to have low weight, as the difference at this bin is not necessarily bad (just because a certain color doesn't
//    * appear in a query histogram, that doesn't mean it is a bad match.) However, differences at bins with small counts in 'idealHistogram'
//    * are important (contribute a lot to the error) because this means that new colors have been introduced, which is bad for "color coherence." */
//  static float WeightedHistogramDifference(const HistogramType& idealHistogram, const HistogramType& queryHistogram);

//  /** Determine how many bins in 'queryHistogram' have non-zero counts while the corresponding bin in 'idealHistogram' has a zero count. The idea is to get a sense of
//    * how many colors are present in the queryHistogram that are not present in the idealHistogram. A low score means that all of the colors in the queryHistogram are
//    * present in the idealHistogram, while a high score means that different colors are present. */
//  static float HistogramCoherence(const HistogramType& idealHistogram, const HistogramType& queryHistogram);

  template <typename TQuadrantHistogram, typename THistogramDifferenceFunctor>
  float QuadrantHistogramDifference(const TQuadrantHistogram& idealHistogram, const TQuadrantHistogram& queryHistogram,
                                    THistogramDifferenceFunctor histogramDifferenceFunctor = [](const typename TQuadrantHistogram::HistogramType& a, const typename TQuadrantHistogram::HistogramType& b)
                                      {
                                      return HistogramDifference(a,b);
                                      })
  {
    float difference = 0.0f;
    for(unsigned int i = 0; i < 4; ++i)
    {
      if(idealHistogram.Properties.Valid[i] != queryHistogram.Properties.Valid[i])
      {
        std::stringstream ss;
        ss << "QuadrantHistograms do not have the same Valid value for quadrant " << i;
        throw std::runtime_error(ss.str());
      }

      if(idealHistogram.Properties.Valid[i] && queryHistogram.Properties.Valid[i])
      {
        difference += histogramDifferenceFunctor(idealHistogram.Histograms[i], queryHistogram.Histograms[i]);
      }
    }

    return difference;
  }

  template <typename THistogram>
  float HistogramCoherence(const THistogram& idealHistogram, const THistogram& queryHistogram)
  {
    if(idealHistogram.size() != queryHistogram.size())
    {
      std::stringstream ss;
      ss << "HistogramCoherence: Histograms must be the same size! idealHistogram is " << idealHistogram.size() << " while queryHistogram is " << queryHistogram.size();
      throw std::runtime_error(ss.str());
    }

    unsigned int numberOfNewBins = 0;
    for(unsigned int bin = 0; bin < idealHistogram.size(); ++bin)
    {
      if((idealHistogram[bin] == 0) && (queryHistogram[bin] > 0))
      {
        numberOfNewBins++;
      }
    }

    return numberOfNewBins;
  }

  template <typename THistogram>
  float WeightedHistogramDifference(const THistogram& idealHistogram, const THistogram& queryHistogram)
  {
    // assert(TBinValue is a signed type)
    if(idealHistogram.size() != queryHistogram.size())
    {
      std::stringstream ss;
      ss << "Histograms must be the same size! idealHistogram is " << idealHistogram.size()
         << " while queryHistogram is " << queryHistogram.size();
      throw std::runtime_error(ss.str());
    }

    float difference = 0.0f;

    for(unsigned int bin = 0; bin < idealHistogram.size(); ++bin)
    {
      // Large values get a small weight
      //float weight = 1.0f/static_cast<float>(idealHistogram[bin]); // cannot do this, because some bins are 0
      float weight = 100.0f - static_cast<float>(idealHistogram[bin]); // this 100 is arbitrary, but if our patches have only a few hundred
          // pixels, and there are ~50 bins, then there should never be a bin with > 100 count.
      difference += weight * fabs(static_cast<float>(idealHistogram[bin]) - static_cast<float>(queryHistogram[bin]));
    }

    return difference;
  }

  template <typename THistogram>
  float HistogramDifference(const THistogram& histogram1, const THistogram& histogram2)
  {
    static_assert(std::is_signed<typename THistogram::value_type>::value, "In HistogramDifference, T must be a signed!");

    // assert(TBinValue is a signed type)
    if(histogram1.size() != histogram2.size())
    {
      std::stringstream ss;
      ss << "Histograms must be the same size! histogram1 is " << histogram1.size() << " while histogram2 is " << histogram2.size();
      throw std::runtime_error(ss.str());
    }

    //float difference = 0.0f;
    typename THistogram::BinValueType difference = 0;
    for(unsigned int bin = 0; bin < histogram1.size(); ++bin)
    {
      //difference += fabs(static_cast<float>(histogram1[bin]) - static_cast<float>(histogram2[bin]));
      difference += fabs(histogram1[bin] - histogram2[bin]);
    }

    return difference;
  }

  /** This function counts how many pixels in 'histogramToCheck' are in bins that aren't occupied in 'idealHistogram' */
  template <typename THistogram>
  unsigned int CountNewColors(const THistogram& idealHistogram, const THistogram& histogramToCheck)
  {
    static_assert(std::is_integral<typename THistogram::value_type>::value, "In CountNewColors, T must be an integral type!");
    if(idealHistogram.size() != histogramToCheck.size())
    {
      std::stringstream ss;
      ss << "Histograms must be the same size! idealHistogram is " << idealHistogram.size()
         << " while histogramToCheck is " << histogramToCheck.size();
      throw std::runtime_error(ss.str());
    }

    unsigned int numberOfPixelsWithNewColor = 0;
    unsigned int numberOfNewColorBins = 0;

    int numberToConsiderEmpty = 3; // If this number of fewer pixels are in a bin, we consider it empty
    // This is just to prevent random outliers from causing a bin to be erroneously empty/occupied.
    for(unsigned int bin = 0; bin < idealHistogram.size(); ++bin)
    {
      if((idealHistogram[bin] <= numberToConsiderEmpty) && (histogramToCheck[bin] > numberToConsiderEmpty))
      {
        numberOfPixelsWithNewColor += (histogramToCheck[bin] - idealHistogram[bin]);
        numberOfNewColorBins++;
      }
    }

//    std::cout << "There are " << numberOfNewColorBins
//              << " new occupied bins (numberOfPixelsWithNewColor = " << numberOfPixelsWithNewColor << ")." << std::endl;

    return numberOfPixelsWithNewColor;
  }

  template <typename THistogram>
  float HistogramIntersection(const THistogram& histogram1, const THistogram& histogram2)
  {
    if(histogram1.size() != histogram2.size())
    {
      std::cerr << "Histograms must be the same size!" << std::endl;
      return 0;
    }

    typedef typename THistogram::BinValueType BinValueType;
    BinValueType totalIntersection = 0.0f;
    for(unsigned int bin = 0; bin < histogram1.size(); ++bin)
    {
      // The casts to float are necessary other wise the integer division always ends up = 0 !
      BinValueType binCount1 = histogram1[bin];
      BinValueType binCount2 = histogram2[bin];
      //std::cout << "frequency1: " << frequency1 << std::endl;
      //std::cout << "frequency2: " << frequency2 << std::endl;
      BinValueType intersection = std::min(binCount1, binCount2);
      //std::cout << "intersection: " << intersection << std::endl;
      totalIntersection += intersection;
    }

    // Why is this of histogram1 instead of histogram2? Are we assuming they have the same total bin sum?
    float totalFrequency = std::accumulate(histogram1.begin(), histogram1.begin() + histogram1.size(), 0.0f); // This 0.0f is the "initial value of the sum". It has to be a float, or the sum will be an integer sum.
    //std::cout << "totalFrequency: " << totalFrequency << std::endl;
    float normalizedIntersection = static_cast<float>(totalIntersection) / totalFrequency;
    //std::cout << "normalizedIntersection: " << normalizedIntersection << std::endl;
    return normalizedIntersection;
  }

}

#endif
