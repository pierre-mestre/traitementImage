#include "tpConnectedComponents.h"
#include <cmath>
#include <algorithm>
#include <tuple>
#include <vector>
#include <map>
#include <stack>
using namespace cv;
using namespace std;

#define WHITE 255


/**
    Performs a labeling of image connected component with 4 connectivity
    with a depth-first exploration.
    Any non zero pixel of the image is considered as present.
*/
cv::Mat ccLabel(cv::Mat image)
{
    Mat res = Mat::zeros(image.rows, image.cols, CV_32SC1); // 32 int image
    /********************************************
                YOUR CODE HERE
    *********************************************/
    stack<cv::Point2i> s;
    Mat visited = Mat::zeros(image.rows, image.cols, CV_8UC1); // 32 int image
    std::vector<cv::Point2i> voisins = {{0,1}, {1,0}, {0,-1}, {-1,0}};
    int ccNumber = 1;

    for (int i = 0; i < res.rows; i++) {
        for (int j = 0; j < res.cols; j++) {
            // Valeur du pixel > 0 et non visité
            if(image.at<int>(i,j) > 0 && !visited.at<uchar>(i,j)) {
                s.push({i,j});
                while(!s.empty()) {
                    cv::Point2i p = s.top();
                    s.pop();

                    for (cv::Point2i v : voisins) {
                        v += p;
                        if(image.at<int>(v.x,v.y) > 0 && !visited.at<uchar>(v.x,v.y)) {
                            visited.at<uchar>(v.x,v.y) = (uchar)1;
                            s.push(v);
                        }
                    }
                    res.at<int>(p.x,p.y) = ccNumber;
                    visited.at<uchar>(p.x,p.y) = (uchar)1;
                }
                ccNumber++;
            }
        }
    }
    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}

/**
    Deletes the connected components (4 connectivity) containg less than size pixels.
*/
cv::Mat ccAreaFilter(cv::Mat image, int size)
{
    Mat res = Mat::zeros(image.rows, image.cols, image.type());
    assert(size>0);


    /********************************************
                YOUR CODE HERE
    *********************************************/
    //
    int ccNumber = 100000;
    int maxCcNumber = 0;

    cv::Mat imageLabel = ccLabel(image);

    // Search max cc
    for (int i = 0; i < res.rows; ++i) {
        for (int j = 0; j < res.cols; ++j) {
            int val = imageLabel.at<int>(i,j);
            if(val > maxCcNumber) {
                maxCcNumber = val;
            }

            if(val < ccNumber && val != 0) {
                ccNumber = val;
            }
        }
    }

    while(ccNumber < maxCcNumber) {
        int sizeCal = 0;
        for (int i = 0; i < res.rows; ++i) {
            for (int j = 0; j < res.cols; ++j) {
                int val = imageLabel.at<int>(i,j);
                if(val == ccNumber) {
                    sizeCal++;
                }
            }
        }

        if(sizeCal > size) {
            for (int i = 0; i < res.rows; i++) {
                for (int j = 0; j < res.cols; j++) {
                    int val = imageLabel.at<int>(i,j);
                    if(val == ccNumber) {
                        res.at<int>(i,j) = image.at<int>(i,j);
                    }
                }
            }
        }
        ccNumber++;
    }

    int count = 0;
    for (int i = 0; i < res.rows; i++) {
        for (int j = 0; j < res.cols; j++) {
            if(res.at<int>(i,j) != 0)
                count++;
        }
    }

    cout << "FINAL NUMBER OF PIXEL: " << count << endl;

    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}


/**
    Performs a labeling of image connected component with 4 connectivity using a
    2 pass algorithm.
    Any non zero pixel of the image is considered as present.
*/

int Find(int x, map<int, int> eq){
    while(eq[x] != x) {
        x = eq[x];
    }
    return x;
}

void Union(int x, int y, std::map<int, int> &eq){
    eq[Find(x, eq)] = Find(y, eq);
}

cv::Mat ccTwoPassLabel(cv::Mat image)
{
    Mat res = Mat::zeros(image.rows, image.cols, CV_32SC1); // 32 int image
    /********************************************
                YOUR CODE HERE
    *********************************************/
    std::vector<Point2i> voisins = {{0, -1}, {-1, 0}};
    int label = 1;
    std::map<int, int> eq;

    // 1ere passe
    for(int i = 0; i < res.rows; ++i) {
        for(int j = 0; j < res.cols; ++j) {
            if (image.at<int>(i, j) == 0) continue;

            Point2i p = {i,j};
            std::vector<Point2i> l;

            for(Point2i v : voisins) {
                v += p;
                if (image.at<int>(v.x, v.y) > 0)
                    l.push_back(v);
            }

            if (!l.empty()) {
                int min = 100000;

                for (Point2i v : l) {
                    if (res.at<int>(v.x, v.y) > 0 && min > res.at<int>(v.x, v.y))
                        min = res.at<int>(v.x, v.y);
                }

                res.at<int>(i, j) = min;

                for (Point2i v  : l) {
                    int px = res.at<int>(v.x, v.y);
                    Union(min, px, eq);
                }
            } else {
                eq[label] = label;
                res.at<int>(i, j) = label;
                label++;
            }
        }
    }
    // 2nde passe
    for(int i = 0; i < res.rows; ++i) {
        for(int j = 0; j < res.cols; ++j) {
            if (image.at<int>(i, j) > 0)
                res.at<int>(i,j) = Find(res.at<int>(i,j), eq);
        }
    }


    /********************************************
                END OF YOUR CODE
    *********************************************/
    return res;
}