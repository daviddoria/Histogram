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
#include <ITKHelpers/ITKHelpers.h>

int main()
{
  //typedef itk::VectorImage<unsigned char, 2> ImageType;
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;

  unsigned int imageWidth = 100;
  
  ImageType::IndexType corner = {{0,0}};
  ImageType::SizeType size = {{imageWidth, imageWidth}};
  ImageType::RegionType region(corner, size);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  // Create a random image
  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

  while(!imageIterator.IsAtEnd())
    {
    ImageType::PixelType pixel(image->GetNumberOfComponentsPerPixel());
    for(unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
    {
      pixel[i] = rand() % 255;
    }

    imageIterator.Set(pixel);
    ++imageIterator;
    }

  TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
  TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;

  unsigned int numberOfBinsPerComponent = 20;
  typedef int BinValueType;

  unsigned int patchRadius = 7;
  //unsigned int patchRadius = 3;

  //itk::Size<2> patchSize = {{patchRadius*2 + 1, patchRadius*2 + 1}};

  float totalDifference = 0;
  
  typedef HistogramGenerator<BinValueType> HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;

  for(unsigned int i = 0; i < 10000; ++i)
  {
    // I'm not sure why the *2's are necessary... but it was creating patches outside of the image without them
    int center1x = rand() % (imageWidth - patchRadius*2 - 1) + patchRadius;
    int center1y = rand() % (imageWidth - patchRadius*2 - 1) + patchRadius;
    itk::Index<2> center1 = {{center1x, center1y}};
    itk::ImageRegion<2> region1 = ITKHelpers::GetRegionInRadiusAroundPixel(center1, patchRadius);

    int center2x = rand() % (imageWidth - patchRadius*2 - 1) + patchRadius;
    int center2y = rand() % (imageWidth - patchRadius*2 - 1) + patchRadius;
    itk::Index<2> center2 = {{center2x, center2y}};
    itk::ImageRegion<2> region2 = ITKHelpers::GetRegionInRadiusAroundPixel(center2, patchRadius);
    
    HistogramType histogram1 = HistogramGeneratorType::ComputeImageHistogram1D(image.GetPointer(),
                                                         region1,
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

    HistogramType histogram2 = HistogramGeneratorType::ComputeImageHistogram1D(image.GetPointer(),
                                                         region2,
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

    float difference = HistogramDifferences::HistogramDifference(histogram1, histogram2);
    totalDifference += difference;
  }

  std::cout << "Difference: " << totalDifference << std::endl;
  
  return 0;
}
