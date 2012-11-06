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

static void TestCompute1DConcatenatedHistogramOfMultiChannelImage();

static void TestScalarHistogram();

static void TestHistogramIntersection();

static void TestHistogramDifference();

static void TestWriteHistogram();

static void TestOutputHistogram();

int main()
{
  TestCompute1DConcatenatedHistogramOfMultiChannelImage();
  TestScalarHistogram();
  TestHistogramDifference();
  TestHistogramIntersection();
  TestWriteHistogram();
  TestOutputHistogram();
  return 0;
}

void TestCompute1DConcatenatedHistogramOfMultiChannelImage()
{
  // Single channel
  {
  typedef itk::Image<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::IndexType corner = {{0,0}};

  ImageType::SizeType size = {{100,100}};

  ImageType::RegionType region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.GetIndex()[0] < 70)
      {
      imageIterator.Set(255);
      }
    else
      {
      imageIterator.Set(0);
      }

    ++imageIterator;
    }

  ImageType::PixelType rangeMin = 0;
  ImageType::PixelType rangeMax = 255;

  unsigned int numberOfBins = 10;
  typedef int BinValueType;
  typedef HistogramGenerator<BinValueType> HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;

  HistogramType histogram = HistogramGeneratorType::ComputeImageHistogram1D(image.GetPointer(),
                                                         image->GetLargestPossibleRegion(),
                                                         numberOfBins, rangeMin, rangeMax);

  histogram.Print();
  std::cout << std::endl;
  }

//   // Multi channel VectorImage
//   {
//   typedef itk::VectorImage<unsigned char, 2> ImageType;
//   ImageType::Pointer image = ImageType::New();
//   ImageType::IndexType corner = {{0,0}};
// 
//   ImageType::SizeType size = {{100,100}};
// 
//   ImageType::RegionType region(corner, size);
// 
//   image->SetRegions(region);
//   image->SetNumberOfComponentsPerPixel(3);
//   image->Allocate();
// 
//   itk::ImageRegionIterator<ImageType> imageIterator(image,region);
// 
//   while(!imageIterator.IsAtEnd())
//     {
//     ImageType::PixelType pixel(image->GetNumberOfComponentsPerPixel());
//     if(imageIterator.GetIndex()[0] < 70)
//       {
//       for(unsigned int i = 0; i < pixel.GetSize(); ++i)
//         {
//         pixel[i] = 255;
//         }
//       }
//     else
//       {
//       for(unsigned int i = 0; i < pixel.GetSize(); ++i)
//         {
//         pixel[i] = 0;
//         }
//       }
//     imageIterator.Set(pixel);
//     ++imageIterator;
//     }
// 
//   TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
//   TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;
// 
//   unsigned int numberOfBinsPerComponent = 10;
//   typedef int BinValueType;
//   Histogram<BinValueType>::HistogramType histogram = Histogram<BinValueType>::Compute1DConcatenatedHistogramOfMultiChannelImage(image.GetPointer(),
//                                                          image->GetLargestPossibleRegion(),
//                                                          numberOfBinsPerComponent, rangeMin, rangeMax);
// 
//   Histogram<BinValueType>::OutputHistogram(histogram);
//   std::cout << std::endl;
//   }

  // Multi channel Image<CovariantVector>
  {
  typedef itk::Image<itk::CovariantVector<unsigned char, 3>, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::IndexType corner = {{0,0}};

  ImageType::SizeType size = {{100,100}};

  ImageType::RegionType region(corner, size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

  while(!imageIterator.IsAtEnd())
    {
    ImageType::PixelType pixel(image->GetNumberOfComponentsPerPixel());
    if(imageIterator.GetIndex()[0] < 70)
      {
      for(unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
        {
        pixel[i] = 255;
        }
      }
    else
      {
      for(unsigned int i = 0; i < pixel.GetNumberOfComponents(); ++i)
        {
        pixel[i] = 0;
        }
      }
    imageIterator.Set(pixel);
    ++imageIterator;
    }

  TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
  TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;

  unsigned int numberOfBinsPerComponent = 10;
  typedef int BinValueType;
  typedef HistogramGenerator<BinValueType> HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;
  HistogramType histogram = HistogramGeneratorType::ComputeImageHistogram1D(image.GetPointer(),
                                                         image->GetLargestPossibleRegion(),
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

  histogram.Print();
  std::cout << std::endl;
  }
}

void TestScalarHistogram()
{
  typedef float ValueType;
  std::vector<ValueType> values;
  values.push_back(1);
  values.push_back(2);
  values.push_back(3);

  unsigned int numberOfBins = 5;
  ValueType rangeMin = 0.0f;
  ValueType rangeMax = 4.0f;

  typedef int BinValueType;
  typedef HistogramGenerator<BinValueType> HistogramGeneratorType;
  typedef HistogramGeneratorType::HistogramType HistogramType;
  HistogramType histogram = HistogramGeneratorType::ScalarHistogram(values, numberOfBins,
                                                            rangeMin, rangeMax);

  histogram.Print();
  std::cout << std::endl;
}

void TestHistogramDifference()
{
  typedef int BinValueType;
  typedef Histogram<BinValueType> HistogramType;

  HistogramType histogram1;
  histogram1.push_back(1);
  histogram1.push_back(2);
  histogram1.push_back(3);

  HistogramType histogram2;
  histogram2.push_back(1);
  histogram2.push_back(2);
  histogram2.push_back(4);
  float difference = HistogramDifferences::HistogramDifference(histogram1, histogram2);

  std::cout << "difference: " << difference << std::endl;
}

void TestHistogramIntersection()
{
  typedef int BinValueType;
  typedef Histogram<BinValueType> HistogramType;
  HistogramType histogram1;
  histogram1.push_back(1);
  histogram1.push_back(2);
  histogram1.push_back(3);

  HistogramType histogram2;
  histogram2.push_back(1);
  histogram2.push_back(2);
  histogram2.push_back(4);
  float intersection = HistogramDifferences::HistogramIntersection(histogram1, histogram2);

  std::cout << "intersection: " << intersection << std::endl;
}

void TestWriteHistogram()
{
  typedef int BinValueType;
  typedef Histogram<BinValueType> HistogramType;
  HistogramType histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(3);

  histogram.Write("histogram.txt");
}

void TestOutputHistogram()
{
  typedef int BinValueType;
  typedef Histogram<BinValueType> HistogramType;
  HistogramType histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(3);

  histogram.Print();
}
