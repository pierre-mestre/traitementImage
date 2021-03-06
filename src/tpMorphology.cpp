#include "tpMorphology.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <limits>
#include "common.h"
using namespace cv;
using namespace std;


/**
    Compute a median filter of the input float image.
    The filter window is a square of (2*size+1)*(2*size+1) pixels.

    Values outside image domain are ignored.

    The median of a list l of n>2 elements is defined as:
     - l[n/2] if n is odd 
     - (l[n/2-1]+l[n/2])/2 is n is even 
*/
Mat median(Mat image, int k)
{
    Mat res = image.clone();
    assert(k>0);
    /********************************************
                YOUR CODE HERE
    *********************************************/

    for(int i = 0; i < image.rows; ++i){
        for (int j = 0; j < image.cols; ++j){
            std::vector<float> m;
            for(int v = i - k; v <= i + k; ++v){
                for (int w = j - k; w <= j + k; ++w){
                    if(v >=0 && v < image.rows && w >= 0 && w < image.cols)
                        m.push_back(image.at<float>(v,w));
                }
            }
            float med(0);
            auto size = m.size();
            std::sort(m.begin(), m.begin()+ size);
            if(m.size()%2)
                med = m[size/2];
            else
                med = (float)(m[size/2-1]+m[size/2])/2;
            res.at<float>(i,j) = med;
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Compute the erosion of the input float image by the given structuring element.
    Pixel outside the image are supposed to have value 1.
*/
Mat erode(Mat image, Mat structuringElement) {
    Mat res = image.clone();
    /********************************************
                YOUR CODE HERE
    *********************************************/
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            res.at<float>(i, j) = 255 - res.at<float>(i, j);
        }
    }
    Mat res_inv = dilate(res, structuringElement);
    for (int i = 0; i < image.rows; ++i) {
        for (int j = 0; j < image.cols; ++j) {
            res_inv.at<float>(i, j) = 255 - res_inv.at<float>(i, j);
        }
    }
    res = res_inv;
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Compute the dilation of the input float image by the given structuring element.
     Pixel outside the image are supposed to have value 0
*/
Mat dilate(Mat image, Mat kernel)
{
    Mat res = Mat::zeros(1,1,CV_32FC1);
    /********************************************
                YOUR CODE HERE
        hint : 1 line of code is enough
    *********************************************/
    res = image.clone();
    int hr = kernel.rows/2;
    int hc = kernel.cols/2;
    for(int i = 0; i < image.rows; ++i){
        for (int j = 0; j < image.cols; ++j){
            std::vector<float> m;
            for(int k = 0; k < kernel.rows; ++k){
                for (int l = 0; l < kernel.cols; ++l){
                    if((i + hr - k) >= 0 && (i + hr - k) < image.rows && (j + hc - l) >= 0 && (j + hc - l) < image.cols) {
                        m.push_back((float) kernel.at<float>(k,l) * image.at<float>(i - k + hr,j - l + hc));
                    }
                }
            }
            res.at<float>(i,j) = (float) *max_element(m.begin(), m.end());
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Compute the opening of the input float image by the given structuring element.
*/
Mat open(Mat image, Mat structuringElement)
{

    Mat res = Mat::zeros(1,1,CV_32FC1);
    /********************************************
                YOUR CODE HERE
        hint : 1 line of code is enough
    *********************************************/
    res = dilate(erode(image, structuringElement),structuringElement);
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Compute the closing of the input float image by the given structuring element.
*/
Mat close(Mat image, Mat structuringElement)
{

    Mat res = Mat::zeros(1,1,CV_32FC1);
    /********************************************
                YOUR CODE HERE
        hint : 1 line of code is enough
    *********************************************/
    res = erode(dilate(image,structuringElement), structuringElement);
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Compute the morphological gradient of the input float image by the given structuring element.
*/
Mat morphologicalGradient(Mat image, Mat structuringElement)
{

    Mat res = Mat::zeros(1,1,CV_32FC1);
    /********************************************
                YOUR CODE HERE
        hint : 1 line of code is enough
    *********************************************/
    res = dilate(image, structuringElement) - erode(image, structuringElement);
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

