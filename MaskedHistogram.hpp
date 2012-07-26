
std::vector<float> Compute1DHistogramOfMultiChannelMaskedImage(const FloatVectorImageType* image, const itk::ImageRegion<2>& imageRegion, const Mask* mask, const itk::ImageRegion<2>& maskRegion, const unsigned int numberOfBins)
{
  // Compute the histogram for each channel separately
  std::vector<HistogramType::Pointer> channelHistograms = ComputeHistogramsOfMaskedRegion(image, imageRegion, mask, maskRegion, numberOfBins);

  unsigned int binsPerHistogram = channelHistograms[0]->GetSize(0);
  std::vector<float> histogram(binsPerHistogram*image->GetNumberOfComponentsPerPixel());

  for(unsigned int channel = 0; channel < image->GetNumberOfComponentsPerPixel(); ++channel)
    {
    for(unsigned int bin = 0; bin < binsPerHistogram; ++bin)
      {
      histogram[channel*binsPerHistogram + bin] = channelHistograms[channel]->GetFrequency(bin);
      }
    }

  return histogram;
}
