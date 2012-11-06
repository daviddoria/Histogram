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

#include "MaskedHistogramGenerator.h"

// ITK
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

// Submodules
#include <Helpers/TypeTraits.h>

static void TestComputeMaskedImage1DHistogram();

int main()
{
  TestComputeMaskedImage1DHistogram();

  return 0;
}

void TestComputeMaskedImage1DHistogram()
{
//  // Single channel
//  {
//  typedef itk::Image<unsigned char, 2> ImageType;
//  ImageType::Pointer image = ImageType::New();
//  ImageType::IndexType corner = {{0,0}};

//  ImageType::SizeType size = {{100,100}};

//  ImageType::RegionType region(corner, size);

//  image->SetRegions(region);
//  image->Allocate();

//  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

//  while(!imageIterator.IsAtEnd())
//    {
//    if(imageIterator.GetIndex()[0] < 70)
//      {
//      imageIterator.Set(255);
//      }
//    else
//      {
//      imageIterator.Set(0);
//      }

//    ++imageIterator;
//    }

//  ImageType::PixelType rangeMin = 0;
//  ImageType::PixelType rangeMax = 255;

//  unsigned int numberOfBins = 10;
//  typedef int BinValueType;
//  typedef Histogram<BinValueType>::HistogramType HistogramType;

//  HistogramType histogram = Histogram<BinValueType>::ComputeImageHistogram1D(image.GetPointer(),
//                                                         image->GetLargestPossibleRegion(),
//                                                         numberOfBins, rangeMin, rangeMax);

//  Histogram<BinValueType>::OutputHistogram(histogram);
//  std::cout << std::endl;
//  }

   // Multi channel VectorImage
   {
   typedef itk::VectorImage<unsigned char, 2> ImageType;
   ImageType::Pointer image = ImageType::New();
   ImageType::IndexType corner = {{0,0}};

   ImageType::SizeType size = {{100,100}};

   ImageType::RegionType region(corner, size);

   image->SetRegions(region);
   image->SetNumberOfComponentsPerPixel(3);
   image->Allocate();

   Mask::Pointer mask = Mask::New();
   mask->SetRegions(region);
   mask->Allocate();

   itk::ImageRegionIterator<ImageType> imageIterator(image,region);

   while(!imageIterator.IsAtEnd())
   {
     ImageType::PixelType pixel(image->GetNumberOfComponentsPerPixel());
     if(imageIterator.GetIndex()[0] < 70)
     {
       for(unsigned int i = 0; i < pixel.GetSize(); ++i)
       {
         pixel[i] = 255;
       }
     }
     else
     {
       for(unsigned int i = 0; i < pixel.GetSize(); ++i)
       {
         pixel[i] = 0;
       }
     }
     imageIterator.Set(pixel);
     ++imageIterator;
   }

//   TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
//   TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;

   ImageType::PixelType rangeMins;
   rangeMins.SetSize(image->GetNumberOfComponentsPerPixel());
   rangeMins.Fill(0);

   ImageType::PixelType rangeMaxs;
   rangeMaxs.SetSize(image->GetNumberOfComponentsPerPixel());
   rangeMaxs.Fill(255);

   unsigned int numberOfBinsPerComponent = 10;
   typedef int BinValueType;

   itk::ImageRegion<2> imageRegion = image->GetLargestPossibleRegion();
   itk::ImageRegion<2> maskRegion = image->GetLargestPossibleRegion();

   typedef MaskedHistogramGenerator<BinValueType> HistogramGeneratorType;
   typedef HistogramGeneratorType::HistogramType HistogramType;

   bool allowOutside = false;
   unsigned char maskValue = mask->GetValidValue();
   HistogramType histogram =
       HistogramGeneratorType::ComputeMaskedImage1DHistogram(image.GetPointer(), imageRegion,
                                                      mask.GetPointer(), maskRegion,
                                                      numberOfBinsPerComponent, rangeMins, rangeMaxs,
                                                      allowOutside, maskValue);

   histogram.Print();
   std::cout << std::endl;
   }

//  // Multi channel Image<CovariantVector>
//  {
//  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;
//  ImageType::Pointer image = ImageType::New();
//  ImageType::IndexType corner = {{0,0}};

//  ImageType::SizeType size = {{100,100}};

//  ImageType::RegionType region(corner, size);

//  image->SetRegions(region);
//  image->Allocate();

//  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

//  while(!imageIterator.IsAtEnd())
//    {
//    ImageType::PixelType pixel(image->GetNumberOfComponentsPerPixel());
//    if(imageIterator.GetIndex()[0] < 70)
//      {
//      for(unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
//        {
//        pixel[i] = 255;
//        }
//      }
//    else
//      {
//      for(unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
//        {
//        pixel[i] = 0;
//        }
//      }
//    imageIterator.Set(pixel);
//    ++imageIterator;
//    }

//  TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
//  TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;

//  unsigned int numberOfBinsPerComponent = 10;
//  typedef int BinValueType;
//  Histogram<BinValueType>::HistogramType histogram = Histogram<BinValueType>::ComputeImageHistogram1D(image.GetPointer(),
//                                                         image->GetLargestPossibleRegion(),
//                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

//  Histogram<BinValueType>::OutputHistogram(histogram);
//  std::cout << std::endl;
//  }
}
