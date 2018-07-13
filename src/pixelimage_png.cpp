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

#include "pixelimage_png.h"

#include <lilyengine/filesystem.h>

#include <png.h>

#include <csetjmp>
#include <cstring>

struct PixelImage_PngReadState
{
    const std::string *data;
    size_t index;
};

static void pixelImage_pngErrorFunc(png_structp png_ptr, const char *message)
{
    longjmp(png_jmpbuf(png_ptr), 1);
}

static void pixelImage_pngReadFunc(png_structp png_ptr, unsigned char *dst, size_t count)
{
    PixelImage_PngReadState *io = (PixelImage_PngReadState*)png_get_io_ptr(png_ptr);
    size_t i = 0;
    while(i < count && io->index < io->data->size()) {
        dst[i++] = (*io->data)[io->index++];
    }
}

ExPop::PixelImage<uint8_t> *pixelImageLoadPNG(const std::string &data)
{
    ExPop::PixelImage_Dimension width = 0;
    ExPop::PixelImage_Dimension height = 0;
    ExPop::PixelImage_Dimension channelCount = 0;
    unsigned char colorType = 0;
    size_t rowBytes = 0;

    ExPop::PixelImage<uint8_t> *newImage = nullptr;

    if(data.size() > 8) {

        if(!png_sig_cmp((const unsigned char *)&data[0], 0, 8)) {

            png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

            if(png_ptr) {

                png_infop info_ptr = png_create_info_struct(png_ptr);

                if(info_ptr) {

                    unsigned char *rawData = nullptr;
                    unsigned char **rowPtrs = nullptr;

                    if(!setjmp(png_jmpbuf(png_ptr))) {

                        struct PixelImage_PngReadState readState;
                        readState.data = &data;
                        readState.index = 0;

                        png_set_read_fn(png_ptr, &readState, pixelImage_pngReadFunc);

                        png_set_error_fn(png_ptr, nullptr, pixelImage_pngErrorFunc, pixelImage_pngErrorFunc);

                        png_read_info(png_ptr, info_ptr);

                        png_set_palette_to_rgb(png_ptr);
                        png_set_tRNS_to_alpha(png_ptr);
                        png_set_expand_gray_1_2_4_to_8(png_ptr);

                        // FIXME: In the future, we can maybe just do
                        // arbitrary bit depths. But right now we'll
                        // just stick to coercing everything to be
                        // 8-bit per channel.
                        png_set_strip_16(png_ptr);

                        png_read_update_info(png_ptr, info_ptr);

                        width = png_get_image_width(png_ptr, info_ptr);
                        height = png_get_image_height(png_ptr, info_ptr);
                        colorType = png_get_color_type(png_ptr, info_ptr);

                        rowBytes = png_get_rowbytes(png_ptr, info_ptr);
                        size_t byteCount = rowBytes * height;
                        rawData = new unsigned char[byteCount];

                        memset(rawData, 0, byteCount);

                        rowPtrs = new unsigned char *[height];

                        for(ExPop::PixelImage_Dimension i = 0; i < height; i++) {
                            rowPtrs[i] = rawData + (rowBytes * i);
                        }

                        png_read_image(png_ptr, rowPtrs);
                    }

                    // Determine channel count from format.
                    if(rawData) {
                        channelCount = 1;
                        if(colorType & PNG_COLOR_MASK_COLOR) {
                            channelCount += 2;
                        }
                        if(colorType & PNG_COLOR_MASK_ALPHA) {
                            channelCount += 1;
                        }
                    }

                    if(channelCount) {
                        newImage = new ExPop::PixelImage<uint8_t>(width, height, channelCount);
                    }

                    // Move data into the actual PixelImage.
                    if(newImage) {
                        for(ExPop::PixelImage_Dimension y = 0; y < height; y++) {
                            for(ExPop::PixelImage_Dimension x = 0; x < width; x++) {
                                for(ExPop::PixelImage_Dimension c = 0; c < channelCount; c++) {
                                    ExPop::PixelImage<uint8_t>::PixelValueType &t = newImage->getData(x, y, c);
                                    t.value = (rawData[(x * channelCount) + (y * rowBytes) + c]);
                                }
                            }
                        }
                    }

                    delete[] rawData;
                    delete[] rowPtrs;

                    png_destroy_info_struct(png_ptr, &info_ptr);
                }

                png_destroy_read_struct(&png_ptr, nullptr, nullptr);
            }
        }
    }

    return newImage;
}

ExPop::PixelImage<uint8_t> *pixelImageLoadPNGFromFile(const std::string &filename)
{
    std::string fileData = ExPop::FileSystem::loadFileString(filename);
    if(fileData.size()) {
        return pixelImageLoadPNG(fileData);
    }
    return nullptr;
}

struct PixelImage_PngWriteState
{
    std::vector<uint8_t> outBuffer;
};

void pixelImagePNGWriteFunc(png_structp png_ptr, unsigned char *src, size_t count)
{
    PixelImage_PngWriteState *io = (PixelImage_PngWriteState*)png_get_io_ptr(png_ptr);
    io->outBuffer.resize(io->outBuffer.size() + count);
    uint8_t *newDataStart = &io->outBuffer[io->outBuffer.size() - count];
    memcpy(newDataStart, src, count);
}

uint8_t *pixelImageSavePNG(
    const ExPop::PixelImage<uint8_t> &img, size_t *length)
{
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    uint8_t *ret = nullptr;
    *length = 0;
    uint8_t *tmpImageData = nullptr;
    uint8_t **rowPtrs = nullptr;

    if(png_ptr) {
        png_infop info_ptr = png_create_info_struct(png_ptr);
        if(info_ptr) {
            if(!setjmp(png_jmpbuf(png_ptr))) {
                png_set_error_fn(png_ptr, nullptr, pixelImage_pngErrorFunc, pixelImage_pngErrorFunc);

                PixelImage_PngWriteState writeState;

                png_set_write_fn(png_ptr, &writeState, pixelImagePNGWriteFunc, nullptr);

                uint32_t colorType = 0;
                switch(img.getChannelCount()) {
                    case 1:
                        colorType = PNG_COLOR_TYPE_GRAY;
                        break;
                    case 2:
                        colorType = PNG_COLOR_TYPE_GRAY_ALPHA;
                        break;
                    case 3:
                        colorType = PNG_COLOR_TYPE_RGB;
                        break;
                    case 4:
                        colorType = PNG_COLOR_TYPE_RGBA;
                        break;
                    default:
                        longjmp(png_jmpbuf(png_ptr), 1);
                        break;
                }

                png_set_IHDR(png_ptr, info_ptr,
                    img.getWidth(), img.getHeight(),
                    8, colorType,
                    PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE,
                    PNG_FILTER_TYPE_BASE);

                png_write_info(png_ptr, info_ptr);

                tmpImageData = new uint8_t[img.getWidth() * img.getHeight() * img.getChannelCount()];
                for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
                    for(ExPop::PixelImage_Dimension x = 0; x < img.getWidth(); x++) {
                        for(ExPop::PixelImage_Dimension c = 0; c < img.getChannelCount(); c++) {
                            tmpImageData[(x * img.getChannelCount()) + (y * img.getWidth() * img.getChannelCount()) + c] =
                                img.getData(x, y, c).value;
                        }
                    }
                }

                rowPtrs = new uint8_t *[img.getHeight()];

                for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
                    rowPtrs[y] = &tmpImageData[(y * img.getWidth() * img.getChannelCount())];
                }

                png_write_image(png_ptr, rowPtrs);
                png_write_end(png_ptr, NULL);

                ret = new uint8_t[writeState.outBuffer.size()];
                memcpy(ret, &writeState.outBuffer[0], writeState.outBuffer.size());
                *length = writeState.outBuffer.size();
            }
            delete[] rowPtrs;
            delete[] tmpImageData;
            png_destroy_info_struct(png_ptr, &info_ptr);
        }
        png_destroy_write_struct(&png_ptr, nullptr);
    }

    return ret;
}

bool pixelImageSavePNGToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename)
{
    size_t imgLen = 0;
    uint8_t *data = pixelImageSavePNG(img, &imgLen);

    bool ret = (0 == ExPop::FileSystem::saveFile(filename, (char*)data, imgLen));

    delete[] data;

    return ret;
}
