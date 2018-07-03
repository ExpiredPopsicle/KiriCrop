#pragma once

#include <lilyengine/pixelimage/pixelimage.h>

ExPop::PixelImage<uint8_t> *pixelImageLoadPNG(const std::string &data);
ExPop::PixelImage<uint8_t> *pixelImageLoadPNGFromFile(const std::string &filename);
uint8_t *pixelImageSavePNG(
    const ExPop::PixelImage<uint8_t> &img, size_t *length);
bool pixelImageSavePNGToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename);


