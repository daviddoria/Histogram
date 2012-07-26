
template <typename TImage>
std::vector<float> Compute1DHistogramOfMultiChannelMaskedImage(const TImage* image,
                                                               const itk::ImageRegion<2>& imageRegion,
                                                               Mask* mask, const itk::ImageRegion<2>& maskRegion,
                                                               const unsigned int numberOfBins);