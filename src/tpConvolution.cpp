
#include "tpConvolution.h"
#include <cmath>
#include <algorithm>
#include <tuple>
using namespace cv;
using namespace std;
/**
    Compute a mean filter of size 2k+1.

    Pixel values outside of the image domain are supposed to have a zero value.
*/
cv::Mat meanFilter(cv::Mat image, int k){
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    float mul = (float)(2 * k + 1);

    for(int i = 0; i < image.rows; ++i){
        for (int j = 0; j < image.cols; ++j){
            float m(0);
            for(int v = i - k; v <= i + k; ++v){
                for (int w = j - k; w <= j + k; ++w){
                    if(v >=0 && v < image.rows && w >= 0 && w < image.cols)
                        m += image.at<float>(v,w);
                }
            }
            res.at<float>(i,j) = m/pow(mul, 2);
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/

    return res;

}

/**
    Compute the convolution of a float image by kernel.
    Result has the same size as image.
    
    Pixel values outside of the image domain are supposed to have a zero value.
*/
Mat convolution(Mat image, cv::Mat kernel)
{
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    int hr = kernel.rows/2;
    int hc = kernel.cols/2;
    for(int i = 0; i < image.rows; ++i){
        for (int j = 0; j < image.cols; ++j){
            float m(0);
            for(int k = 0; k < kernel.rows; ++k){
                for (int l = 0; l < kernel.cols; ++l){
                    if((i + hr - k) >= 0 && (i + hr - k) < image.rows && (j + hc - l) >= 0 && (j + hc - l) < image.cols) {
                        m += (float) kernel.at<float>(k,l) * image.at<float>(i - k + hr,j - l + hc);
                    }
                }
            }
            res.at<float>(i,j) = m;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/

    return res;
}

/**
    Compute the sum of absolute partial derivative according to Sobel's method
*/
cv::Mat edgeSobel(cv::Mat image)
{
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    auto t = image.type();
    float arr_v1[3] = {1,2,1};
    float arr_v2[3] = {-1,0,1};
    float arr_h1[3] = {-1,0,1};
    float arr_h2[3] = {1,2,1};

    res = abs(convolution(convolution(image, cv::Mat(3, 1, t, arr_v1)), cv::Mat(1, 3, t, arr_h1))) +
            abs(convolution(convolution(image, cv::Mat(3, 1, t, arr_v2)), cv::Mat(1, 3, t, arr_h2)));
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Value of a centered gaussian of variance (scale) sigma at point x.
*/
float gaussian(float x, float sigma2)
{
    return 1.0/(2*M_PI*sigma2)*exp(-x*x/(2*sigma2));
}

/**
    Performs a bilateral filter with the given spatial smoothing kernel 
    and a intensity smoothing of scale sigma_r.

*/
cv::Mat bilateralFilter(cv::Mat image, cv::Mat kernel, float sigma_r)
{
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/

    int hr = kernel.rows/2;
    int hc = kernel.cols/2;
    for(int i = 0; i < image.rows; ++i){
        for (int j = 0; j < image.cols; ++j){
            float m(0), n(0);
            for(int k = 0; k < kernel.rows; ++k){
                for (int l = 0; l < kernel.cols; ++l){
                    if((i + hr - k) >= 0 && (i + hr - k) < image.rows && (j + hc - l) >= 0 && (j + hc - l) < image.cols) {
                        float o = gaussian(image.at<float>(i + hr - k, j + hc - l) - image.at<float>(i,j), sigma_r*sigma_r);
                        m += (float) kernel.at<float>(k,l) * image.at<float>(i - k + hr,j - l + hc) * o;
                        n += (float) kernel.at<float>(k,l) * o;
                    }
                }
            }
            res.at<float>(i,j) = m/n;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}
