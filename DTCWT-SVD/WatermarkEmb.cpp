#include <iostream>
#include "WatermarkEmb.h"

#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\imgproc\imgproc.hpp>

using namespace std;
using namespace cv;

int main() {
	Mat originalImage = imread("E:/DTCWT-SVD/DTCWT-SVD/src/001.jpg");
	int watermark = 1;
	char* saveImagePath = (char*) "E:/DTCWT-SVD/DTCWT-SVD/src/out.tiff";
	int T = 200;
	int a = 3;
	int lamda = 3;

	Mat out = watermarkEmbImg(originalImage, watermark, saveImagePath, T, a, lamda);
	imwrite("E:/DTCWT-SVD/DTCWT-SVD/src/out.jpg", out);

	int level = 3;

	int waterGet = watermarkGetImg(out, level);
	cout << "The watermark is: " << waterGet;

	system("pause");
	return 0;
}

Mat watermarkEmbImg(Mat original_image, int watermark, char *saveImagePath, int threshold, int a, int lamda) {
	
	int original_image_height = original_image.rows;
	int original_image_width = original_image.cols;

	//获取图像的YUV数据
	BYTE *R_arr = new BYTE[original_image_height*original_image_width]();
	BYTE *G_arr = new BYTE[original_image_height*original_image_width]();
	BYTE *B_arr = new BYTE[original_image_height*original_image_width]();
	BYTE *Y_arr = new BYTE[original_image_height*original_image_width]();
	BYTE *U_arr = new BYTE[original_image_height*original_image_width]();
	BYTE *V_arr = new BYTE[original_image_height*original_image_width]();
	double T[3][3] = { { 0.2126,0.7152,0.0722 },{ -0.1146,-0.3854,0.5000 },{ 0.5000,-0.4542,-0.0468 } };    //RGB转YUV参数
																										//获取图像RGB值
	int m = 0;
	for (int i = 0; i < original_image_height; i++)
	{
		for (int j = 0; j < original_image_width; j++)
		{
			int Blue = original_image.at<Vec3b>(i, j)[0];
			int Green = original_image.at<Vec3b>(i, j)[1];
			int Red = original_image.at<Vec3b>(i, j)[2];

			R_arr[m] = Red;
			G_arr[m] = Green;
			B_arr[m] = Blue;
			m = m + 1;
		}
	}

	//RGB转成YUV
	for (int x = 0; x < original_image_height*original_image_width; x++)
	{
		Y_arr[x] = std::round(T[0][0] * R_arr[x] + T[0][1] * G_arr[x] + T[0][2] * B_arr[x]);
		U_arr[x] = std::round(T[1][0] * R_arr[x] + T[1][1] * G_arr[x] + T[1][2] * B_arr[x]) + 127;
		V_arr[x] = std::round(T[2][0] * R_arr[x] + T[2][1] * G_arr[x] + T[2][2] * B_arr[x]) + 127;
	}

	//嵌入
	int *Out_Uarr = new int[original_image_height*original_image_width]();
	imageTools* tools = new imageTools(watermark,threshold, a, lamda, 0, original_image_height, original_image_width);
	tools->AddWaterMark(U_arr, Out_Uarr);
	tools->WaterEncodeDeinit();

	//输出图像
	int *OutR_arr = new int[original_image_height*original_image_width]();
	int *OutG_arr = new int[original_image_height*original_image_width]();
	int *OutB_arr = new int[original_image_height*original_image_width]();
	double T1[3][3] = { { 1,0,1.4075 },{ 1,-0.3455,-0.7169 },{ 1,1.779,0 } };

	//YUV转RGB
	double num_min = 0;
	double num_max = 255;
	for (int x = 0; x < original_image_height*original_image_width; x++)
	{
		double r = T1[0][0] * Y_arr[x] + T1[0][1] * (Out_Uarr[x] - 127) + T1[0][2] * (V_arr[x] - 127);
		double g = T1[1][0] * Y_arr[x] + T1[1][1] * (Out_Uarr[x] - 127) + T1[1][2] * (V_arr[x] - 127);
		double b = T1[2][0] * Y_arr[x] + T1[2][1] * (Out_Uarr[x] - 127) + T1[2][2] * (V_arr[x] - 127);
		if (r > num_max)
			r = 255;
		if (r < num_min)
			r = 0;
		if (g > num_max)
			g = 255;
		if (g < num_min)
			g = 0;
		if (b > num_max)
			b = 255;
		if (b < num_min)
			b = 0;
		int r1 = (int)std::round(r);
		int g1 = (int)std::round(g);
		int b1 = (int)std::round(b);
		OutR_arr[x] = r1;
		OutG_arr[x] = g1;
		OutB_arr[x] = b1;
	}

	//生成图像    
	Mat embImage = cv::Mat(original_image_height, original_image_width, CV_8UC3);
	int desb = 0;
	int desg = 0;
	int desr = 0;

	int desm = 0;
	for (int i = 0; i < original_image_height; i++)
	{

		for (int j = 0; j < original_image_width; j++)
		{

			desb = OutB_arr[desm];
			desg = OutG_arr[desm];
			desr = OutR_arr[desm];
			embImage.at<Vec3b>(i, j)[0] = desb;
			embImage.at<Vec3b>(i, j)[1] = desg;
			embImage.at<Vec3b>(i, j)[2] = desr;
			desm = desm + 1;
		}
	}
	//cv::imwrite(SaveImagePath, embImage);

	//delete util;
	delete[]R_arr;
	delete[]G_arr;
	delete[]B_arr;
	delete[]Y_arr;
	delete[]U_arr;
	delete[]V_arr;
	delete[]Out_Uarr;
	delete tools;
	delete[]OutR_arr;
	delete[]OutG_arr;
	delete[]OutB_arr;

	return embImage;
}

int watermarkGetImg(Mat markedImage, int level) {
	return 1;
}


