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

#include "HistogramGenerator.h"
#include "HistogramDifferences.hpp"

// ITK
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

// Submodules
#include <Helpers/TypeTraits.h>
#include <ITKHelpers/itkNormImageAdaptor.h>

int main()
{
  typedef itk::Image<itk::CovariantVector<float, 2>, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::IndexType corner = {{0,0}};
  ImageType::SizeType size = {{100,100}};
  ImageType::RegionType region(corner, size);
  image->SetRegions(region);
  image->Allocate();

  typedef ImageType::PixelType::RealValueType ScalarType;

  typedef itk::NormImageAdaptor<ImageType, ScalarType>
      NormImageAdaptorType;
  typename NormImageAdaptorType::Pointer imageAdaptor = NormImageAdaptorType::New();
  imageAdaptor->SetImage(image);

  ScalarType rangeMin = 0;
  ScalarType rangeMax = 255;

  unsigned int numberOfBins = 10;
  typedef int BinValueType;
  typedef HistogramGenerator<BinValueType> HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;

  HistogramType histogram =
      HistogramGeneratorType::ComputeImageHistogram1D(imageAdaptor.GetPointer(),
                                                      imageAdaptor->GetLargestPossibleRegion(),
                                                      numberOfBins, rangeMin, rangeMax);

  histogram.Print();
  std::cout << std::endl;

  return 0;
}
