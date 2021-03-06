// ---------------------------------------------------------------------------
//
//   Kiri's Cropping Tool
//
//   Copyright (c) 2018 Kiri Jolly
//     http://expiredpopsicle.com
//     expiredpopsicle@gmail.com
//
// ---------------------------------------------------------------------------
//
//   This software is provided 'as-is', without any express or implied
//   warranty. In no event will the authors be held liable for any
//   damages arising from the use of this software.
//
//   Permission is granted to anyone to use this software for any
//   purpose, including commercial applications, and to alter it and
//   redistribute it freely, subject to the following restrictions:
//
//   1. The origin of this software must not be misrepresented; you must
//      not claim that you wrote the original software. If you use this
//      software in a product, an acknowledgment in the product
//      documentation would be appreciated but is not required.
//
//   2. Altered source versions must be plainly marked as such, and must
//      not be misrepresented as being the original software.
//
//   3. This notice may not be removed or altered from any source
//      distribution.
//
//
// -------------------------- END HEADER -------------------------------------

#pragma once

#include <lilyengine/pixelimage/pixelimage.h>

inline uint8_t *pixelImageSaveJPEG(
    const ExPop::PixelImage<uint8_t> &img, size_t *length);
bool pixelImageSaveJPEGToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename);
ExPop::PixelImage<uint8_t> *pixelImageLoadJPEG(const std::string &data);
ExPop::PixelImage<uint8_t> *pixelImageLoadJPEGFromFile(const std::string &filename);
