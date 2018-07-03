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

#include "pixelimage_jpeg.h"

#include <lilyengine/filesystem.h>

#include <jpeglib.h>

#include <csetjmp>
#include <cstring>

inline uint8_t *pixelImageSaveJPEG(
    const ExPop::PixelImage<uint8_t> &img, size_t *length)
{
    jpeg_compress_struct compressStruct = {0};

    jpeg_error_mgr errorManager;
    jpeg_std_error(&errorManager);
    compressStruct.err = &errorManager;

    jpeg_create_compress(&compressStruct);
    uint8_t *outBuffer = nullptr;
    unsigned long outSize = 0;
    jpeg_mem_dest(&compressStruct, &outBuffer, &outSize);

    compressStruct.image_width = img.getWidth();
    compressStruct.image_height = img.getHeight();
    compressStruct.input_components = img.getChannelCount();

    if(img.getChannelCount() == 3) {
        compressStruct.in_color_space = JCS_RGB;
    } else if(img.getChannelCount() == 1) {
        compressStruct.in_color_space = JCS_GRAYSCALE;
    } else {
        // Whatever this is, it's not something we can save here.
        jpeg_destroy_compress(&compressStruct);
        return nullptr;
    }

    jpeg_set_defaults(&compressStruct);
    jpeg_set_quality(&compressStruct, 99, true);
    jpeg_start_compress(&compressStruct, true);

    for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
        uint8_t *rowPtr = (uint8_t*)&img.getData(0, y, 0);
        jpeg_write_scanlines(
            &compressStruct,
            &rowPtr,
            1);
    }

    jpeg_finish_compress(&compressStruct);
    jpeg_destroy_compress(&compressStruct);

    *length = outSize;
    return outBuffer;
}

bool pixelImageSaveJPEGToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename)
{
    size_t imgLen = 0;
    uint8_t *data = pixelImageSaveJPEG(img, &imgLen);

    bool ret = (0 == ExPop::FileSystem::saveFile(filename, (char*)data, imgLen));

    free(data);

    return ret;
}

struct PixelImage_JPEGErrorHandler
{
    jpeg_error_mgr libjpegErrorManager; // Must be first.
    jmp_buf returnJump;
};

void pixelImage_jpegError(j_common_ptr decompressStruct)
{
    PixelImage_JPEGErrorHandler *realErrorHandler =
        (PixelImage_JPEGErrorHandler*)decompressStruct->err;

    longjmp(realErrorHandler->returnJump, 1);
}

ExPop::PixelImage<uint8_t> *pixelImageLoadJPEG(const std::string &data)
{
    jpeg_decompress_struct decompressStruct = {0};

    PixelImage_JPEGErrorHandler errorManager = {0};
    jpeg_std_error(&errorManager.libjpegErrorManager);
    decompressStruct.err = &errorManager.libjpegErrorManager;
    errorManager.libjpegErrorManager.error_exit = pixelImage_jpegError;

    ExPop::PixelImage<uint8_t> *ret = nullptr;

    if(!setjmp(errorManager.returnJump)) {

        jpeg_create_decompress(&decompressStruct);
        jpeg_mem_src(&decompressStruct, (uint8_t*)&data[0], data.size());
        jpeg_read_header(&decompressStruct, true);
        jpeg_start_decompress(&decompressStruct);

        ret = new ExPop::PixelImage<uint8_t>(
            decompressStruct.output_width,
            decompressStruct.output_height,
            decompressStruct.output_components);

        while(decompressStruct.output_scanline < decompressStruct.output_height) {

            uint8_t *scanlinePtr =
                (uint8_t*)&ret->getData(0, decompressStruct.output_scanline, 0);

            jpeg_read_scanlines(
                &decompressStruct,
                &scanlinePtr, 1);
        }

        jpeg_finish_decompress(&decompressStruct);
    }

    jpeg_destroy_decompress(&decompressStruct);

    return ret;
}

ExPop::PixelImage<uint8_t> *pixelImageLoadJPEGFromFile(const std::string &filename)
{
    std::string fileData = ExPop::FileSystem::loadFileString(filename);
    if(fileData.size()) {
        return pixelImageLoadJPEG(fileData);
    }
    return nullptr;
}
