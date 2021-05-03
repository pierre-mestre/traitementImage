#include "tpGeometry.h"
#include <cmath>
#include <algorithm>
#include <tuple>
using namespace cv;
using namespace std;

/**
    Transpose the input image,
    ie. performs a planar symmetry according to the
    first diagonal (upper left to lower right corner).
*/
Mat transpose(Mat image)
{
    Mat res = Mat::zeros(image.cols, image.rows, CV_32FC1);
    /********************************************
                YOUR CODE HERE
    hint: consider a non square image
    *********************************************/
    for(int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            res.at<int>(i, j) = image.at<int>(j, i);
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Compute the value of a nearest neighbour interpolation
    in image Mat at position (x,y)
*/
float interpolate_nearest(Mat image, float y, float x)
{
    float v=0;
    /********************************************
                YOUR CODE HERE
    *********************************************/
    v = image.at<float>((int)floor(x+0.5),(int)floor(y+0.5));
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return v;
}


/**
    Compute the value of a bilinear interpolation in image Mat at position (x,y)
*/
float interpolate_bilinear(Mat image, float y, float x)
{
    float v=0;
    /********************************************
                YOUR CODE HERE
    *********************************************/
    int xi = floor(x);
    int yi = floor(y);

    float n1 = image.at<float>(xi+1,yi+1);
    float n2 = image.at<float>(xi,yi+1);
    float n3 = image.at<float>(xi+1,yi);
    float n4 = image.at<float>(xi,yi);

    float a = (float)(x-xi)/(float)((xi+1)-xi);
    float b = (float)(y-yi)/(float)((yi+1)-yi);
    v = (float)((float)a*b*n1 + (float)(1-a)*b*n2 + (float)a*(1-b)*n3 + (float)(1-a)*(1-b)*n4);
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return v;
}

/**
    Multiply the image resolution by a given factor using the given interpolation method.
    If the input size is (h,w) the output size shall be ((h-1)*factor, (w-1)*factor)
*/
Mat expand(Mat image, int factor, float(* interpolationFunction)(cv::Mat image, float y, float x))
{
    assert(factor>0);
    Mat res = Mat::zeros((image.rows-1)*factor,(image.cols-1)*factor,CV_32FC1);

    /********************************************
                YOUR CODE HERE
    *********************************************/
    int rows = image.rows-1;
    int cols = image.cols-1;

    for(int i = 0; i < rows*factor; ++i){
        for (int j = 0; j < cols*factor; ++j){
            res.at<float>(i,j) =  interpolationFunction(image,(float)j/factor,(float)i/factor);
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Performs a rotation of the input image with the given angle (clockwise) and the given interpolation method.
    The center of rotation is the center of the image.

    Ouput size depends of the input image size and the rotation angle.

    Output pixels that map outside the input image are set to 0.
*/
Mat rotate(Mat image, float angle, float(* interpolationFunction)(cv::Mat image, float y, float x))
{
    Mat res = Mat::zeros(1,1,CV_32FC1);
    /********************************************
                YOUR CODE HERE
    hint: to determine the size of the output, take
    the bounding box of the rotated corners of the 
    input image.
    *********************************************/
    
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;

}