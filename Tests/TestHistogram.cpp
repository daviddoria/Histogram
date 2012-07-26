#include "Histogram.h"

// ITK
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkVectorImage.h"

// Submodules
#include "Helpers/TypeTraits.h"

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
  Histogram::HistogramType histogram = Histogram::Compute1DConcatenatedHistogramOfMultiChannelImage(image.GetPointer(),
                                                         image->GetLargestPossibleRegion(),
                                                         numberOfBins, rangeMin, rangeMax);

  Histogram::OutputHistogram(histogram);
  std::cout << std::endl;
  }

  // Multi channel
  {
  typedef itk::VectorImage<unsigned char, 2> ImageType;
  ImageType::Pointer image = ImageType::New();
  ImageType::IndexType corner = {{0,0}};

  ImageType::SizeType size = {{100,100}};

  ImageType::RegionType region(corner, size);

  image->SetRegions(region);
  image->SetNumberOfComponentsPerPixel(3);
  image->Allocate();

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

  TypeTraits<ImageType::PixelType>::ComponentType rangeMin = 0;
  TypeTraits<ImageType::PixelType>::ComponentType rangeMax = 255;

  unsigned int numberOfBinsPerComponent = 10;
  Histogram::HistogramType histogram = Histogram::Compute1DConcatenatedHistogramOfMultiChannelImage(image.GetPointer(),
                                                         image->GetLargestPossibleRegion(),
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

  Histogram::OutputHistogram(histogram);
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
  std::vector<float> histogram = Histogram::ScalarHistogram(values, numberOfBins,
                                                            rangeMin, rangeMax);

  Histogram::OutputHistogram(histogram);
  std::cout << std::endl;
}

void TestHistogramDifference()
{
  Histogram::HistogramType histogram1;
  histogram1.push_back(1);
  histogram1.push_back(2);
  histogram1.push_back(3);

  Histogram::HistogramType histogram2;
  histogram2.push_back(1);
  histogram2.push_back(2);
  histogram2.push_back(4);
  float difference = Histogram::HistogramDifference(histogram1, histogram2);

  std::cout << "difference: " << difference << std::endl;
}

void TestHistogramIntersection()
{
  Histogram::HistogramType histogram1;
  histogram1.push_back(1);
  histogram1.push_back(2);
  histogram1.push_back(3);

  Histogram::HistogramType histogram2;
  histogram2.push_back(1);
  histogram2.push_back(2);
  histogram2.push_back(4);
  float intersection = Histogram::HistogramIntersection(histogram1, histogram2);

  std::cout << "intersection: " << intersection << std::endl;
}

void TestWriteHistogram()
{
  Histogram::HistogramType histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(3);

  Histogram::WriteHistogram(histogram, "histogram.txt");
}

void TestOutputHistogram()
{
  Histogram::HistogramType histogram;
  histogram.push_back(1);
  histogram.push_back(2);
  histogram.push_back(3);

  Histogram::OutputHistogram(histogram);
}
