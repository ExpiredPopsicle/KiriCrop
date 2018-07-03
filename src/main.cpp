#include <lilyengine/utils.h>

#include <png.h>
#include <jpeglib.h>

#include <iostream>

#include <setjmp.h>

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

                        png_set_read_fn(png_ptr, &readState, readFunc);

                        png_set_error_fn(png_ptr, nullptr, pixelImage_pngErrorFunc, pixelImage_pngErrorFunc);

                        png_read_info(png_ptr, info_ptr);


                        // png_set_gray_to_rgb(png_ptr);
                        png_set_palette_to_rgb(png_ptr);
                        png_set_tRNS_to_alpha(png_ptr);

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

void writeFunc(png_structp png_ptr, unsigned char *src, size_t count)
{
    PixelImage_PngWriteState *io = (PixelImage_PngWriteState*)png_get_io_ptr(png_ptr);
    // size_t i = 0;
    // while(i < count && io->index < io->data->size()) {
    //     dst[i++] = (*io->data)[io->index++];
    // }
    io->outBuffer.resize(io->outBuffer.size() + count);
    uint8_t *newDataStart = &io->outBuffer[io->outBuffer.size() - count];
    memcpy(newDataStart, src, count);
}

inline uint8_t *pixelImageSavePNG(
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

                png_set_write_fn(png_ptr, &writeState, writeFunc, nullptr);

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

        // size_t bufferSize =
        //     decompressStruct.output_width *
        //     decompressStruct.output_height *
        //     decompressStruct.output_components;

        // uint8_t *buffer = new uint8_t[bufferSize];

        ret = new ExPop::PixelImage<uint8_t>(
            decompressStruct.output_width,
            decompressStruct.output_height,
            decompressStruct.output_components);

        // memset(buffer, 0, bufferSize);

        while(decompressStruct.output_scanline < decompressStruct.output_height) {
            // uint8_t *scanlinePtr = &buffer[
            //     decompressStruct.output_width * decompressStruct.output_scanline];
            uint8_t *scanlinePtr =
                (uint8_t*)&ret->getData(0, decompressStruct.output_scanline, 0);

            jpeg_read_scanlines(
                &decompressStruct,
                &scanlinePtr, 1);
        }

        jpeg_finish_decompress(&decompressStruct);
    }

    jpeg_destroy_decompress(&decompressStruct);

    // delete[] buffer;

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

ExPop::PixelImage<uint8_t> *pixelImageLoad(const std::string &fileData)
{
    ExPop::PixelImage<uint8_t> *ret = nullptr;
    if(fileData.size()) {
        ret = pixelImageLoadJPEG(fileData);
        if(!ret)
            ret = pixelImageLoadPNG(fileData);
        if(!ret)
            ret = ExPop::pixelImageLoadTGA(&fileData[0], fileData.size());
    }
    return ret;
}

ExPop::PixelImage<uint8_t> *pixelImageLoadFromFile(const std::string &filename)
{
    std::string fileData = ExPop::FileSystem::loadFileString(filename);
    ExPop::PixelImage<uint8_t> *ret = nullptr;
    if(fileData.size()) {
        return pixelImageLoad(fileData);
    }
    return ret;
}

bool pixelImageSaveToFile(
    const ExPop::PixelImage<uint8_t> &img, const std::string &filename)
{
    bool ret = false;
    if(ExPop::stringEndsWith<char>(".png", filename)) {
        ret = pixelImageSavePNGToFile(img, filename);
    } else if(ExPop::stringEndsWith<char>(".tga", filename)) {
        ret = pixelImageSaveTGAToFile(img, filename);
    } else if(ExPop::stringEndsWith<char>(".jpg", filename) ||
        ExPop::stringEndsWith<char>(".jpeg", filename))
    {
        ret = pixelImageSaveJPEGToFile(img, filename);
    }
    return ret;
}

ExPop::PixelImage<uint8_t> *extractLuminance(const ExPop::PixelImage<uint8_t> &img)
{
    ExPop::PixelImage<uint8_t> *ret = new ExPop::PixelImage<uint8_t>(
        img.getWidth(), img.getHeight(), 1);

    if(img.getChannelCount() == 3 || img.getChannelCount() == 4) {

        for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
            for(ExPop::PixelImage_Dimension x = 0; x < img.getWidth(); x++) {

                const ExPop::PixelImage<uint8_t>::PixelValueType &rv =
                    img.getData(x, y, 0);
                const ExPop::PixelImage<uint8_t>::PixelValueType &gv =
                    img.getData(x, y, 1);
                const ExPop::PixelImage<uint8_t>::PixelValueType &bv =
                    img.getData(x, y, 2);
                float rf = rv.getScaledValue<float>();
                float gf = gv.getScaledValue<float>();
                float bf = bv.getScaledValue<float>();

                float lum =
                    0.2126f * rf +
                    0.7152f * gf +
                    0.0722f * bf;

                ret->getData(x, y, 0).setScaledValue<float>(lum);
            }
        }

    } else if(img.getChannelCount() == 1 || img.getChannelCount() == 2) {

        // Just copy the grey channel right over.
        for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
            for(ExPop::PixelImage_Dimension x = 0; x < img.getWidth(); x++) {
                ret->getData(x, y, 0) = img.getData(x, y, 0);
            }
        }

    } else {

        delete ret;
        ret = nullptr;

    }

    return ret;
}

ExPop::PixelImage<uint8_t> *sobelFilter(const ExPop::PixelImage<uint8_t> &img)
{
    ExPop::PixelImage<uint8_t> *ret = new ExPop::PixelImage<uint8_t>(
        img.getWidth(), img.getHeight(), 1);

    assert(img.getChannelCount() == 1);

    const float kernelH[3][3] = {
        {  1,  0, -1 },
        {  2,  0, -2 },
        {  1,  0, -1 }
    };

    const float kernelV[3][3] = {
        {  1,  2,  1 },
        {  0,  0,  0 },
        { -1, -2, -1 }
    };

    for(ExPop::PixelImage_Dimension y = 0; y < img.getHeight(); y++) {
        for(ExPop::PixelImage_Dimension x = 0; x < img.getWidth(); x++) {

            float magH = 0.0f;
            float magV = 0.0f;

            for(int k = 0; k < 3; k++) {
                for(int j = 0; j < 3; j++) {
                    float v = img.getData(
                        x + k - 1,
                        y + j - 1,
                        0,
                        ExPop::PixelImage_EdgeMode_Clamp).getScaledValue<float>();

                    magH += kernelH[k][j] * v;
                    magV += kernelV[k][j] * v;
                }
            }

            float mag = sqrt(magH * magH + magV * magV);

            ret->getData(x, y, 0).setScaledValue(mag);
        }
    }

    return ret;
}

float getRowOrColValue(
    const ExPop::PixelImage<uint8_t> &img,
    int bigAxis,
    ExPop::PixelImage_Dimension axisPos)
{
    ExPop::PixelImage_Dimension axes[2] = {
        img.getWidth(),
        img.getHeight()
    };

    float mag = 0.0f;

    for(ExPop::PixelImage_Dimension pos = 0; pos < axes[!bigAxis]; pos++) {
        ExPop::PixelImage_Dimension fullPos[2];
        fullPos[!bigAxis] = pos;
        fullPos[bigAxis]  = axisPos;

        mag += img.getData(fullPos[0], fullPos[1], 0).getScaledValue<float>();
    }

    return mag;
}

int main(int argc, char *argv[])
{
    ExPop::CommandlineParser cmdParser(argv[0]);

    // Dimensions.
    ExPop::PixelImage_Dimension width = 512;
    cmdParser.addVariableHandler("width", &width);
    cmdParser.setParameterAlias("width", "w");

    ExPop::PixelImage_Dimension height = 512;
    cmdParser.addVariableHandler("height", &height);
    cmdParser.setParameterAlias("height", "h");

    // Output filename.
    std::string outputFilename;
    cmdParser.addVariableHandler("out", &outputFilename);
    cmdParser.setParameterAlias("out", "o");

    // Input filename handler.
    std::vector<std::string> inputFilenames;
    cmdParser.addHandler<std::string>(
        "",
        [&inputFilenames](const std::string &fname){
            inputFilenames.push_back(fname);
        });

    cmdParser.setDoc(
        "Kiri's Cropping Tool 0.1",
        "[options] <imagefile> [<imagefile>...]",
        R"(This tool attemps to resize and crop an image to the most complex
region to match the aspect ratio of the given size. The region
complexity is considered to be the region with the most edges as
detected by a Sobel edge detection filter done against a luminance map
of the input images for RGB and RGBA inputs, or against the image
itself for Grayscale and Grayscale + alpha images.)",
        R"(Copyright (c) 2018 Kiri Jolly)");

    cmdParser.setParameterDoc("width", "Output image width. Defaults to 512.", "integer");
    cmdParser.setParameterDoc("height", "Output image height. Defaults to 512.", "integer");
    cmdParser.setParameterDoc("out",
        R"(Output image filename. Defaults to the same filename as the input.
        Output name may not be specified separately when multiple
        files are specified for processing.)", "name");

    if(!cmdParser.handleCommandline(argc, argv)) {
        return cmdParser.getErrorFlag();
    }

    if(inputFilenames.size() == 0) {
        std::cerr << cmdParser.getHelpText() << std::endl;
        return 1;
    }

    if(outputFilename.size() != 0 && inputFilenames.size() != 1) {
        std::cerr << "Error: Multiple input files specified when using --out." << std::endl;
        return 1;
    }

    if(height <= 0 || width <= 0) {
        std::cerr << "Error: Invalid size specified." << std::endl;
        return 1;
    }

    float desiredAspectRatio = float(width) / float(height);

    for(size_t i = 0; i < inputFilenames.size(); i++) {

        ExPop::PixelImage<uint8_t> *img = pixelImageLoadFromFile(inputFilenames[i]);

        std::cout << "Processing: " << inputFilenames[i] << std::endl;

        if(img) {
            pixelImageSaveToFile(*img, "out.tga");
        } else {
            std::cerr << "Error: Failed to load image: " << inputFilenames[i] << std::endl;
            std::cerr << "       Unsupported image format or type." << std::endl;
            return 1;
        }


        ExPop::PixelImage<uint8_t> *lum = extractLuminance(*img);
        ExPop::PixelImage<uint8_t> *sobel = sobelFilter(*lum);

        pixelImageSaveToFile(*lum, "lum.png");
        pixelImageSaveToFile(*sobel, "sobel.png");

        float inputAspectRatio =
            float(img->getWidth()) /
            float(img->getHeight());

        std::cout << "Input aspect ratio: " << inputAspectRatio << std::endl;
        std::cout << "Desired aspect ratio: " << desiredAspectRatio << std::endl;

        ExPop::PixelImage_Dimension newHeight = 0;
        ExPop::PixelImage_Dimension newWidth = 0;
        int bigAxis = 0;

        if(inputAspectRatio > desiredAspectRatio) {
            // Width extends beyond desired size.
            newHeight = height;
            newWidth = height * (inputAspectRatio);
            bigAxis = 0;
            std::cout << "Image size now (1): " << newWidth << ", " << newHeight << std::endl;
        } else {
            // Height extends beyond desired size or whatever.
            newHeight = width * (1.0f / inputAspectRatio);
            newWidth = width;
            bigAxis = 1;
            std::cout << "Image size now (2): " << newWidth << ", " << newHeight << std::endl;
        }

        ExPop::PixelImage<uint8_t> *scaledSobel = ExPop::pixelImageScale<uint8_t>(*sobel, newWidth, newHeight);
        ExPop::PixelImage<uint8_t> *scaledImg = ExPop::pixelImageScale<uint8_t>(*img, newWidth, newHeight);


        pixelImageSaveToFile(*scaledImg, "scaledimg.png");
        pixelImageSaveToFile(*scaledSobel, "scaledsobel.png");



        ExPop::PixelImage_Dimension axes[2] = {
            scaledSobel->getWidth(),
            scaledSobel->getHeight()
        };

        ExPop::PixelImage_Dimension newAxes[2] = {
            width,
            height
        };


        float biggestWindowValue = 0.0f;
        ExPop::PixelImage_Dimension biggestWindowStart = 0;
        float currentWindowValue = 0.0f;
        std::vector<float> rowValues;

        for(ExPop::PixelImage_Dimension pos = 0; pos < axes[bigAxis]; pos++) {
            rowValues.push_back(getRowOrColValue(*scaledSobel, bigAxis, pos));
        }

        // Build up initial window.
        for(ExPop::PixelImage_Dimension i = 0; i < newAxes[bigAxis]; i++) {
            currentWindowValue += rowValues[i];
        }

        std::cout << "Initial window: " << currentWindowValue << std::endl;
        biggestWindowStart = 0;
        biggestWindowValue = currentWindowValue;

        std::cout << "newAxes[bigAxis]: " << newAxes[bigAxis] << std::endl;
        std::cout << "axes[bigAxis]:    " << axes[bigAxis] << std::endl;

        for(ExPop::PixelImage_Dimension i = newAxes[bigAxis] + 1; i < axes[bigAxis]; i++) {
            ExPop::PixelImage_Dimension startRow = i - newAxes[bigAxis];
            currentWindowValue -= rowValues[startRow - 1];
            currentWindowValue += rowValues[i];

            if(biggestWindowValue < currentWindowValue) {
                biggestWindowValue = currentWindowValue;
                biggestWindowStart = startRow;
            }
            std::cout << currentWindowValue << std::endl;
        }

        std::cout << "Startrow: " << biggestWindowStart << std::endl;

        ExPop::PixelImage<uint8_t> outputImage(width, height, scaledImg->getChannelCount());
        for(ExPop::PixelImage_Dimension n = 0; n < newAxes[!bigAxis]; n++) {
            for(ExPop::PixelImage_Dimension i = 0; i < newAxes[bigAxis]; i++) {

                ExPop::PixelImage_Dimension dstpos[2];
                dstpos[bigAxis] = i;
                dstpos[!bigAxis] = n;

                ExPop::PixelImage_Dimension srcpos[2];
                srcpos[bigAxis] = i + biggestWindowStart;
                srcpos[!bigAxis] = n;

                for(ExPop::PixelImage_Dimension c = 0; c < scaledImg->getChannelCount(); c++) {
                    ExPop::PixelImage<uint8_t>::PixelValueType &vOut = outputImage.getData(dstpos[0], dstpos[1], c);
                    ExPop::PixelImage<uint8_t>::PixelValueType &vIn = scaledImg->getData(srcpos[0], srcpos[1], c);

                    vOut = vIn;
                }
            }
        }

        if(outputFilename.size() != 0) {
            pixelImageSaveToFile(outputImage, outputFilename);
        } else {
            pixelImageSaveToFile(outputImage, inputFilenames[i]);
        }

        delete sobel;
        delete scaledImg;
        delete scaledSobel;
        delete lum;
        delete img;
    }

    return 0;
}
