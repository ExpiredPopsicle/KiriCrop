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
#include "pixelimage_jpeg.h"

#include <lilyengine/utils.h>

#include <iostream>
#include <csetjmp>

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
        R"(This tool attempts to resize and crop an image to the most complex
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

        ExPop::PixelImage_Dimension newHeight = 0;
        ExPop::PixelImage_Dimension newWidth = 0;
        int bigAxis = 0;

        // Determine which axis we need to crop (the "big axis").
        if(inputAspectRatio > desiredAspectRatio) {

            // Width extends beyond desired size.
            newHeight = height;
            newWidth = height * (inputAspectRatio);
            bigAxis = 0;

        } else {

            // Height extends beyond desired size or whatever.
            newHeight = width * (1.0f / inputAspectRatio);
            newWidth = width;
            bigAxis = 1;

        }

        // Scale the sobel-filtered version of the image and the color image down.
        ExPop::PixelImage<uint8_t> *scaledSobel = ExPop::pixelImageScale<uint8_t>(*sobel, newWidth, newHeight);
        ExPop::PixelImage<uint8_t> *scaledImg = ExPop::pixelImageScale<uint8_t>(*img, newWidth, newHeight);

        // Find the "most interesting" part of the image.
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

        // Figure out per-row complexity values.
        for(ExPop::PixelImage_Dimension pos = 0; pos < axes[bigAxis]; pos++) {
            rowValues.push_back(getRowOrColValue(*scaledSobel, bigAxis, pos));
        }

        // Figure out the complexity value for the initial window.
        for(ExPop::PixelImage_Dimension i = 0; i < newAxes[bigAxis]; i++) {
            currentWindowValue += rowValues[i];
        }

        biggestWindowStart = 0;
        biggestWindowValue = currentWindowValue;

        // Search for the best window.
        for(ExPop::PixelImage_Dimension i = newAxes[bigAxis] + 1; i < axes[bigAxis]; i++) {

            ExPop::PixelImage_Dimension startRow = i - newAxes[bigAxis];
            currentWindowValue -= rowValues[startRow - 1];
            currentWindowValue += rowValues[i];

            if(biggestWindowValue < currentWindowValue) {
                biggestWindowValue = currentWindowValue;
                biggestWindowStart = startRow;
            }
        }

        // Copy the best window to the output image.
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

        // Save the file.
        if(outputFilename.size() != 0) {
            pixelImageSaveToFile(outputImage, outputFilename);
        } else {
            pixelImageSaveToFile(outputImage, inputFilenames[i]);
        }

        // Cleanup.
        delete sobel;
        delete scaledImg;
        delete scaledSobel;
        delete lum;
        delete img;
    }

    return 0;
}
