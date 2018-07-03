#pragma once

#include <lilyengine/pixelimage/pixelimage.h>

inline uint8_t *pixelImageSaveJPEG(
    const ExPop::PixelImage<uint8_t> &img, size_t *length);
bool pixelImageSaveJPEGToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename);
ExPop::PixelImage<uint8_t> *pixelImageLoadJPEG(const std::string &data);
ExPop::PixelImage<uint8_t> *pixelImageLoadJPEGFromFile(const std::string &filename);
