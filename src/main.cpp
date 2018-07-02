#include <lilyengine/utils.h>

#include <png.h>

#include <iostream>

struct PixelImage_PngReadState
{
    const std::string *data;
    size_t index;
};

void pixelImage_pngErrorFunc(png_structp png_ptr, const char *message)
{
    longjmp(png_jmpbuf(png_ptr), 1);
}

void readFunc(png_structp png_ptr, unsigned char *dst, size_t count)
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
    unsigned char bitDepth = 0;
    size_t rowBytes = 0;

    ExPop::PixelImage<uint8_t> *newImage = nullptr;

    if(data.size() > 8) {

        if(!png_sig_cmp((const unsigned char *)&data[0], 0, 8)) {

            png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

            if(png_ptr) {

                png_infop info_ptr = png_create_info_struct(png_ptr);

                if(info_ptr) {

                    std::cout << "READING A PNG" << std::endl;

                    unsigned char *rawData = nullptr;
                    unsigned char **rowPtrs = nullptr;

                    if(!setjmp(png_jmpbuf(png_ptr))) {

                        struct PixelImage_PngReadState readState;
                        readState.data = &data;
                        readState.index = 0;

                        png_set_read_fn(png_ptr, &readState, readFunc);

                        png_set_error_fn(png_ptr, nullptr, pixelImage_pngErrorFunc, pixelImage_pngErrorFunc);

                        png_read_info(png_ptr, info_ptr);


                        png_set_gray_to_rgb(png_ptr);
                        png_set_palette_to_rgb(png_ptr);
                        png_set_tRNS_to_alpha(png_ptr);

                        png_read_update_info(png_ptr, info_ptr);


                        width = png_get_image_width(png_ptr, info_ptr);
                        height = png_get_image_height(png_ptr, info_ptr);
                        colorType = png_get_color_type(png_ptr, info_ptr);
                        bitDepth = png_get_bit_depth(png_ptr, info_ptr);


                        rowBytes = png_get_rowbytes(png_ptr, info_ptr);
                        size_t byteCount = rowBytes * height;
                        rawData = new unsigned char[byteCount];

                        memset(rawData, 0, byteCount);

                        rowPtrs = new unsigned char *[height];

                        for(ExPop::PixelImage_Dimension i = 0; i < height; i++) {
                            rowPtrs[i] = rawData + (rowBytes * i);
                        }

                        // switch(colorType) {
                        //     case PNG_COLOR_TYPE_GRAY:
                        //     case PNG_COLOR_TYPE_GRAY_ALPHA:
                        //         std::cout << "Converting grey to RGB" << std::endl;
                        //         // png_set_gray_to_rgb(png_ptr);
                        //         break;
                        //     case PNG_COLOR_TYPE_PALETTE:
                        //         std::cout << "Converting palette to RGB" << std::endl;
                        //         // png_set_palette_to_rgb(png_ptr);
                        //         break;
                        // }

                        // std::cout << "COLOR TYPE: " << (int)colorType << std::endl;
                        // png_set_tRNS_to_alpha(png_ptr);

                        png_read_image(png_ptr, rowPtrs);
                    }

                    if(rawData) {
                        if(colorType == PNG_COLOR_TYPE_RGB) {
                            channelCount = 3;
                        } else if(colorType == PNG_COLOR_TYPE_RGB_ALPHA) {
                            channelCount = 4;
                        } else {
                            // Can't decode this.
                        }
                    }

                    if(channelCount) {
                        newImage = new ExPop::PixelImage<uint8_t>(width, height, channelCount);
                    }

                    if(newImage) {
                        for(ExPop::PixelImage_Dimension x = 0; x < width; x++) {
                            for(ExPop::PixelImage_Dimension y = 0; y < height; y++) {
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

int main(int argc, char *argv[])
{
    ExPop::PixelImage<uint8_t> *img = pixelImageLoadPNGFromFile(argv[1]);

    std::cout << "Image: " << img << std::endl;

    if(img) {
        pixelImageSaveTGAToFile(*img, "out.tga");
    }

    delete img;

    return 0;
}
