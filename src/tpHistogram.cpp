#include "tpHistogram.h"
#include <cmath>
#include <algorithm>
#include <tuple>
using namespace cv;
using namespace std;

/**
    Inverse a grayscale image with float values.
    for all pixel p: res(p) = 1.0 - image(p)
*/
Mat inverse(Mat image)
{
    // clone original image
    Mat res = image.clone();
      /********************************************
                YOUR CODE HERE
    *********************************************/
    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            res.at<float>(i,j) = 1.0f - image.at<float>(i,j);
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Thresholds a grayscale image with float values.
    for all pixel p: res(p) =
        | 0 if image(p) <= lowT
        | image(p) if lowT < image(p) <= hightT
        | 1 otherwise
*/
Mat threshold(Mat image, float lowT, float highT)
{
    Mat res = image.clone();
    assert(lowT <= highT);
    /********************************************
                YOUR CODE HERE
    *********************************************/
    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            float val = image.at<float>(i,j);
            res.at<float>(i,j) = val < lowT ? 0 : val > highT ? 1.0f : val;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Quantize the input float image in [0,1] in numberOfLevels different gray levels.
    eg. for numberOfLevels = 3 the result should be for all pixel p: res(p) =
        | 0 if image(p) < 1/3
        | 1/2 if 1/3 <= image(p) < 2/3
        | 1 otherwise
*/
Mat quantize(Mat image, int numberOfLevels)
{
    Mat res = image.clone();
    assert(numberOfLevels>0);
    /********************************************
                YOUR CODE HERE
    *********************************************/
    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            float val = image.at<float>(i,j);
            for (int k = 0; k < numberOfLevels; ++k) {
                if((float)k/(float)numberOfLevels <= (float)val && (float)val < (float)(k+1)/(float)numberOfLevels)
                    res.at<float>(i,j) = (float)k/(float)(numberOfLevels-1);
            }
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Normalize a grayscale image with float values
    Target range is [minValue, maxValue].
*/
Mat normalize(Mat image, float minValue, float maxValue)
{
    Mat res = image.clone();
    assert(minValue <= maxValue);
    /********************************************
                YOUR CODE HERE
    *********************************************/
    double min, max;
    minMaxLoc(image, &min, &max);

    float diffAct = max - min;
    float diffRes = maxValue - minValue;

    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            float val = image.at<float>(i,j);
            res.at<float>(i,j) = (val - min) * diffRes / diffAct + minValue;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}



/**
    Equalize image histogram with unsigned char values ([0;255])
*/
Mat equalize(Mat image)
{
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    int pixels = res.rows * res.cols;
    int nbPixels[256] = {0};
    int histo[256] = {0};

    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            nbPixels[(int)res.at<uchar>(i,j)]++;
        }
    }

    histo[0] = nbPixels[0];
    for (int i = 1; i < 256; ++i) {
        histo[i] = histo[i-1] + nbPixels[i];
    }

    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            int val = (int)res.at<uchar>(i,j);
            res.at<uchar>(i,j) = (uchar)((255.0/pixels) * histo[val] + 0.5);
        }
    }


    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;

}

/**
    Compute a binarization of the input float image using an automatic Otsu threshold.
    Input image is of type unsigned char ([0;255])
*/
Mat thresholdOtsu(Mat image)
{
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    int pixels = res.rows * res.cols;
    float sum = 0.0f;
    float sum2 = 0.0f;
    float a = 0.0f;
    float max = 0.0f;
    float thr = 0.0f;
    float histo[256] = {0};

    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            uchar val = (uchar)res.at<uchar>(i,j);
            histo[val]++;
            thr++;
        }
    }

    for (int i = 0; i < 256; ++i) {
        sum += histo[i]*i;
    }

    for (int i = 0; i < 256; ++i) {
        a += histo[i];

        float a2 = pixels - a;
        if(a2 == 0) break; //avoid divide by 0

        sum2 += i * histo[i];

        float moy1 = sum2/a;
        float moy2 = (sum - sum2)/a2;

        float var = a*a2*((moy2-moy1)*(moy2-moy1));

        if(var > max) {
            thr = i;
            max = var;
        }

    }

    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            res.at<uchar>(i,j) = res.at<uchar>(i,j) > thr ? 255 : 0;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}
