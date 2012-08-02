#include "Histogram.h"

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
  
  for(unsigned int i = 0; i < 10000; ++i)
  {
    // I'm not sure why the *2's are necessary... but it was creating patches outside of the image without them
    itk::Index<2> center1 = {{rand() % (imageWidth - patchRadius*2 - 1) + patchRadius,
                             rand() % (imageWidth - patchRadius*2 - 1) + patchRadius}};
    itk::ImageRegion<2> region1 = ITKHelpers::GetRegionInRadiusAroundPixel(center1, patchRadius);

    itk::Index<2> center2 = {{rand() % (imageWidth - patchRadius*2 - 1) + patchRadius,
                             rand() % (imageWidth - patchRadius*2 - 1) + patchRadius}};
    itk::ImageRegion<2> region2 = ITKHelpers::GetRegionInRadiusAroundPixel(center2, patchRadius);
    
    Histogram<BinValueType>::HistogramType histogram1 = Histogram<BinValueType>::ComputeImageHistogram1D(image.GetPointer(),
                                                         region1,
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

    Histogram<BinValueType>::HistogramType histogram2 = Histogram<BinValueType>::ComputeImageHistogram1D(image.GetPointer(),
                                                         region2,
                                                         numberOfBinsPerComponent, rangeMin, rangeMax);

    float difference = Histogram<BinValueType>::HistogramDifference(histogram1, histogram2);
    totalDifference += difference;
  }

  std::cout << "Difference: " << totalDifference << std::endl;
  
  return 0;
}
