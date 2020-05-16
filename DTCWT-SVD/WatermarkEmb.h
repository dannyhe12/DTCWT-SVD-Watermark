#pragma once
#include "imageTools.h"
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\imgproc\imgproc.hpp>

using namespace cv;
using namespace std;

Mat watermarkEmbImg(Mat originalImage, int watermark, char *saveImagePath, int threshold, int a, int lamda);
int watermarkGetImg(Mat markedImage, int level);