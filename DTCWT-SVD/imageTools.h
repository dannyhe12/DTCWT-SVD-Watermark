#pragma once

#include <cstring>
#include <math.h>
#include <vector>
#include <iostream>
#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\opencv.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <iomanip>
#include <csvd.h>

using namespace std;
using namespace splab;
using namespace cv;

#define sqrt2	0.707106781186548	//根号2分之1
typedef unsigned char BYTE;
typedef vector<double> Array;

class imageTools {
public:
	imageTools();
	~imageTools();
	imageTools(int watermark, int threshold,int a,int lamda, int level, int h, int w);

	double m_intensity;
	int m_orgwidth, m_orgheight;                                        //原图像大小
	int m_wwidth, m_wheight;                                        //水印大小

																	//水印图像参数
	double *wh;

	//图像参数
	double *Uh1;     //高频部分
	double *L1L1;

	double *Uh2;
	double *L2L2;

	double *Uh3;
	double *L3L3;

	double *Uhabs;


	//输出
	double *Z1;
	double *Z2;
	double *Z3;

	double *newU2;



	//数组大小参数  用于U通道值修改后 返回原图像时数组的大小
	int mgai_Hhalf2;            //270  这是U第三层修改后此时的大小
	int mgai_Whalf2;            //480
	int mnewL2L2_H;             //540    U第三层往第二层返回时，应有的大小
	int mnewL2L2_W;            //960

	int mgai_Hhalf;            //540    原本L2L2的高和宽
	int mgai_Whalf;            //960


	int mgai_H;                //原始L1L1大小
	int mgai_W;

	int mnewL1L1_H;            //这个是判断L1L1得到后判断能否被4整除后的新L1L1大小
	int mnewL1L1_W;

	int sg0o;
	int sg1o;
	int sh0o;
	int sh1o;


	void colfilter(double X[], int XHeight, int XWidth, double filter[], int filtersize, double out[], bool addflag);
	void coldfilt(double X[], int H, int W, double filter1[], double filter2[], double out[]);
	void colifilt(double X[], int H, int W, double filter1[], double filter2[], double out[]);
	void dtwavexfm2(double in[], int H, int W, double out[], double LoLo[], int level);
	void dtwaveifm2(double in[], int H, int W, double Zin[], double Zout[], int Level);
	void AddWaterMark(BYTE *InU, int *OutU);
	void q2c(double in[], int H, int W, double out[], int  num1, int num2);
	void c2q(double in[], int H, int W, double out[], int  num1, int num2);
	void add2Dtcwt(int watermark, int H, int W, double Uh[]);
	void canbe4(vector<vector<double> > arr, vector<vector<double> > *out);
	void extendH(vector<vector<double> > *arr, int H);
	void extendW(vector<vector<double> > *arr);
	void getNewRow(vector<vector<double> > in, vector<vector<double> > *out);
	void getNewCol(vector<vector<double> > in, vector<vector<double> > *out);
	void getNewRow2(vector<vector<double> > in, vector<vector<double> > *out);    //获取第一行到倒数第二行
	void getNewCol2(vector<vector<double> > in, vector<vector<double> > *out);    //获取第一列到倒数第二列
	int WaterEncodeDeinit();
	vector<Matrix<complex<double>>> rearrange(double Uh[], int H, int W);
	double* arrange(vector<Matrix<complex<double>>> res, int H, int W);
	double getMatrixMean(Matrix<double> in);
	
	
	//取矩阵的子矩阵
	template <typename T>
	T getSubMatrix(T in, int startRow, int endRow, int startCol, int endCol) {
		int row = endRow - startRow + 1;
		int col = endCol - startCol + 1;

		T res(row, col);
		int m = 0, n = 0;
		for (int i = startRow; i <= endRow; i++) {
			n = 0;
			for (int j = startCol; j <= endCol; j++) {
				res[m][n] = in[i][j];
				n++;
			}
			m++;
		}
		return res;
	}

	//将子矩阵按坐标位置放入矩阵
	template <typename S>
	S setMatrixVal(S M, S subM, int startRow, int endRow, int startCol, int endCol) {
		S res = M;
		int m = 0, n = 0;
		for (int i = startRow; i <= endRow; i++) {
			n = 0;
			for (int j = startCol; j <= endCol; j++) {
				res[i][j] = subM[m][n];
				n++;
			}
			m++;
		}
		return res;
	}
};