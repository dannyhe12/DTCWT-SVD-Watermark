#define BOUNDS_CHECK

#include "imageTools.h"

using namespace std;
using namespace cv;

double g0o[19] = { 0.0000706263950892857, 0.0,                 -0.00134190150669643, -0.00188337053571429, 0.00715680803571429,
0.0238560267857143,   -0.0556431361607143,  -0.0516880580357143,   0.299757603236607,
0.559430803571429,     0.299757603236607,   -0.0516880580357143,  -0.0556431361607143,
0.0238560267857143,    0.00715680803571429, -0.00188337053571429, -0.00134190150669643,
0.0, 0.0000706263950892857 };

double g1o[13] = { -0.00175781250000000, 0,                  0.0222656250000000, 0.0468750000000000, -0.0482421875000000,
-0.296875000000000,   0.555468750000000, -0.296875000000000, -0.0482421875000000,
0.0468750000000000,  0.0222656250000000, 0,                 -0.00175781250000000 };

double h0o[13] = { -0.00175781250000000, 0,                  0.0222656250000000, -0.0468750000000000, -0.0482421875000000,
0.296875000000000,   0.555468750000000,  0.296875000000000,  -0.0482421875000000,
-0.0468750000000000,  0.0222656250000000, 0,                  -0.00175781250000000 };

double h1o[19] = { -0.0000706263950892857, 0,                    0.00134190150669643, -0.00188337053571429, -0.00715680803571429,
0.0238560267857143,    0.0556431361607143,  -0.0516880580357143,  -0.299757603236607,
0.559430803571429,    -0.299757603236607,   -0.0516880580357143,   0.0556431361607143,
0.0238560267857143,   -0.00715680803571429, -0.00188337053571429,  0.00134190150669643,
0, -0.0000706263950892857 };
//qshift_b
double g0a[14] = { -0.00455689562847549, -0.00543947593727412, 0.0170252238815540, 0.0238253847949203,
-0.106711804686665,    0.0118660920337970,  0.568810420712123,  0.756145643892523,
0.275295384668882,   -0.117203887699115,  -0.0388728012688278, 0.0346603468448535,
-0.00388321199915849,  0.00325314276365318 };

double g0b[14] = { 0.00325314276365318, -0.00388321199915849, 0.0346603468448535, -0.0388728012688278,
-0.117203887699115,    0.275295384668882,   0.756145643892523,   0.568810420712123,
0.0118660920337970,  -0.106711804686665,   0.0238253847949203,  0.0170252238815540,
-0.00543947593727412, -0.00455689562847549 };

double g1a[14] = { -0.00325314276365318, -0.00388321199915849, -0.0346603468448535, -0.0388728012688278,
0.117203887699115,    0.275295384668882,   -0.756145643892523,   0.568810420712123,
-0.0118660920337970,  -0.106711804686665,   -0.0238253847949203,  0.0170252238815540,
0.00543947593727412, -0.00455689562847549 };

double g1b[14] = { -0.00455689562847549,  0.00543947593727412, 0.0170252238815540, -0.0238253847949203,
-0.106711804686665,   -0.0118660920337970,  0.568810420712123,  -0.756145643892523,
0.275295384668882,    0.117203887699115,  -0.0388728012688278, -0.0346603468448535,
-0.00388321199915849, -0.00325314276365318 };

double h0a[14] = { 0.00325314276365318, -0.00388321199915849, 0.0346603468448535, -0.0388728012688278,
-0.117203887699115,    0.275295384668882,   0.756145643892523,   0.568810420712123,
0.0118660920337970,  -0.106711804686665,   0.0238253847949203,  0.0170252238815540,
-0.00543947593727412, -0.00455689562847549 };

double h0b[14] = { -0.00455689562847549, -0.00543947593727412, 0.0170252238815540, 0.0238253847949203,
-0.106711804686665,    0.0118660920337970,  0.568810420712123,  0.756145643892523,
0.275295384668882,   -0.117203887699115,  -0.0388728012688278, 0.0346603468448535,
-0.00388321199915849,  0.00325314276365318 };

double h1a[14] = { -0.00455689562847549,  0.00543947593727412, 0.0170252238815540, -0.0238253847949203,
-0.106711804686665,   -0.0118660920337970,  0.568810420712123,  -0.756145643892523,
0.275295384668882,    0.117203887699115,  -0.0388728012688278, -0.0346603468448535,
-0.00388321199915849, -0.00325314276365318 };

double h1b[14] = { -0.00325314276365318, -0.00388321199915849, -0.0346603468448535, -0.0388728012688278,
0.117203887699115,    0.275295384668882,   -0.756145643892523,   0.568810420712123,
-0.0118660920337970,  -0.106711804686665,   -0.0238253847949203,  0.0170252238815540,
0.00543947593727412, -0.00455689562847549 };

int watermark, T, a, lamda, getLevel;

imageTools::imageTools(){
}

imageTools::~imageTools(){
}

imageTools::imageTools(int water, int threshold,int a1,int lamda1, int level, int h, int w) {
	
	watermark = water;
	T = threshold;
	a = a1;
	lamda = lamda1;
	getLevel = level;
	m_orgwidth = w;
	m_orgheight = h;

	sg0o = 19;
	sg1o = 13;
	sh0o = 13;
	sh1o = 19;
}

void imageTools::colfilter(double X[], int XHeight, int XWidth, double filter[], int filtersize, double out[], bool addflag)
{
	int size2 = filtersize / 2;
	int length = XHeight + 2 * size2;

	int xe_len = XHeight + 2 * size2;
	int* xe = new int[xe_len];

	for (int i = 0; i < XHeight; i++)
		xe[size2 + i] = i + 1;

	for (int j = 0; j < size2; j++)		//小于min的序号 负数
		xe[size2 - 1 - j] = xe[size2 + j];

	for (int j = 0; j < size2; j++)		//超过max的序号
		xe[size2 + XHeight + j] = xe[size2 + XHeight - 1 - j];

	double* conv = new double[length * XWidth];
	for (int i = 0; i < length; i++)
	{
		int index = xe[i] - 1;
		for (int j = 0; j < XWidth; j++)
			conv[i * XWidth + j] = X[index * XWidth + j];
	}

	int Hout = xe_len - filtersize + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < XWidth; j++)
		{
			double con = 0;
			for (int m = 0; m < filtersize; m++)
				con += conv[(i + m) * XWidth + j] * filter[filtersize - 1 - m];

			if (addflag)						//转置
				out[j * Hout + i] += con;
			else
				out[j * Hout + i] = con;
		}
	}
	delete xe;
	delete conv;
}


void imageTools::coldfilt(double X[], int H, int W, double filter1[], double filter2[], double out[])//默认滤波器长度为14
{
	int filtersize = 14;
	int filtersize2 = filtersize / 2;
	int xe_len = H + 2 * filtersize;
	int *xe = new int[xe_len];

	for (int i = 0; i < H; i++)
		xe[filtersize + i] = i + 1;

	for (int j = 0; j < filtersize; j++)		//小于min的序号 负数
		xe[filtersize - 1 - j] = xe[filtersize + j];

	for (int j = 0; j < filtersize; j++)		//超过max的序号
		xe[filtersize + H + j] = xe[filtersize + H - 1 - j];

	double *filter1o = new double[filtersize2];
	double *filter1e = new double[filtersize2];
	double *filter2o = new double[filtersize2];
	double *filter2e = new double[filtersize2];

	for (int i = 0; i < filtersize2; i++)
	{
		filter1o[i] = filter1[i * 2];
		filter1e[i] = filter1[i * 2 + 1];
		filter2o[i] = filter2[i * 2];
		filter2e[i] = filter2[i * 2 + 1];
	}

	int length = (H + 2 * filtersize) / 4 - 1;
	int *t = new int[length];

	for (int i = 0; i < length; i++)
	{
		t[i] = i * 4 + 6;
	}

	double* conv1 = new double[length * W];
	double* conv2 = new double[length * W];
	double* conv3 = new double[length * W];
	double* conv4 = new double[length * W];

	for (int i = 0; i < length; i++)
	{
		int index1 = t[i] - 1 - 1;
		int index2 = t[i] - 3 - 1;
		int index3 = t[i] - 1;
		int index4 = t[i] - 2 - 1;
		index1 = xe[index1] - 1;
		index2 = xe[index2] - 1;
		index3 = xe[index3] - 1;
		index4 = xe[index4] - 1;
		for (int j = 0; j < W; j++)
		{
			conv1[i*W + j] = X[index1*W + j];
			conv2[i*W + j] = X[index2*W + j];
			conv3[i*W + j] = X[index3*W + j];
			conv4[i*W + j] = X[index4*W + j];
		}
	}

	double sum = 0;
	for (int i = 0; i < filtersize; i++)
		sum += filter1[i] * filter2[i];

	double con1, con2, con3, con4;
	int Hout = length - filtersize2 + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < W; j++)
		{
			con1 = 0;
			con2 = 0;
			con3 = 0;
			con4 = 0;
			for (int m = 0; m < filtersize2; m++)
			{
				con1 += conv1[(i + m) * W + j] * filter1o[filtersize2 - 1 - m];
				con2 += conv2[(i + m) * W + j] * filter1e[filtersize2 - 1 - m];
				con3 += conv3[(i + m) * W + j] * filter2o[filtersize2 - 1 - m];
				con4 += conv4[(i + m) * W + j] * filter2e[filtersize2 - 1 - m];
			}
			if (sum > 0)
			{
				out[j * Hout * 2 + 2 * i] = con1 + con2;		//转置!!!!!这里判断
				out[j * Hout * 2 + 2 * i + 1] = con3 + con4;
			}
			else
			{
				out[j * Hout * 2 + 2 * i + 1] = con1 + con2;	//转置!!!!!这里判断
				out[j * Hout * 2 + 2 * i] = con3 + con4;
			}
		}
	}

	delete filter1o;
	delete filter1e;
	delete filter2o;
	delete filter2e;

	delete conv1;
	delete conv2;
	delete conv3;
	delete conv4;

	delete xe;
	delete t;
}


void imageTools::colifilt(double X[], int H, int W, double filter1[], double filter2[], double out[])
{
	int filtersize = 14;
	int filtersize2 = filtersize / 2;
	int xe_len = H + 2 * filtersize2;
	int* xe = new int[xe_len];

	double* filter1o = new double[filtersize2];
	double* filter1e = new double[filtersize2];
	double* filter2o = new double[filtersize2];
	double* filter2e = new double[filtersize2];

	int length = (H + filtersize) / 2 - 1;
	int* t = new int[length];
	double* conv1 = new double[length * W];
	double* conv2 = new double[length * W];

	for (int i = 0; i < H; i++)
		xe[filtersize2 + i] = i + 1;
	for (int j = 0; j < filtersize2; j++)		//小于min的序号
		xe[filtersize2 - 1 - j] = xe[filtersize2 + j];
	for (int j = 0; j < filtersize2; j++)		//超过max的序号
		xe[filtersize2 + H + j] = xe[filtersize2 + H - 1 - j];

	for (int i = 0; i < filtersize2; i++)
	{
		filter1o[i] = filter1[i * 2];
		filter1e[i] = filter1[i * 2 + 1];
		filter2o[i] = filter2[i * 2];
		filter2e[i] = filter2[i * 2 + 1];
	}

	for (int i = 0; i < length; i++)
		t[i] = i * 2 + 3;

	double sum = 0;
	for (int i = 0; i < filtersize; i++)
		sum += filter1[i] * filter2[i];


	int index1, index2;
	for (int i = 0; i < length; i++)
	{
		if (sum > 0)	//ta tv
		{
			index1 = t[i] - 1;
			index2 = t[i] - 2;
		}
		else
		{
			index1 = t[i] - 2;
			index2 = t[i] - 1;
		}
		index1 = xe[index1] - 1;
		index2 = xe[index2] - 1;

		for (int j = 0; j < W; j++)
		{
			conv1[i * W + j] = X[index2 * W + j];
			conv2[i * W + j] = X[index1 * W + j];
		}
	}

	double con1, con2, con3, con4;
	int Hout = length - filtersize2 + 1;
	for (int i = 0; i < Hout; i++)
	{
		for (int j = 0; j < W; j++)
		{
			con1 = 0;
			con2 = 0;
			con3 = 0;
			con4 = 0;
			for (int m = 0; m < filtersize2; m++)
			{
				con1 += conv1[(i + m) * W + j] * filter1o[filtersize2 - 1 - m];
				con2 += conv2[(i + m) * W + j] * filter2o[filtersize2 - 1 - m];
				con3 += conv1[(i + m) * W + j] * filter1e[filtersize2 - 1 - m];
				con4 += conv2[(i + m) * W + j] * filter2e[filtersize2 - 1 - m];
			}
			out[j * Hout * 4 + 4 * i] += con1;		//转置
			out[j * Hout * 4 + 4 * i + 1] += con2;
			out[j * Hout * 4 + 4 * i + 2] += con3;
			out[j * Hout * 4 + 4 * i + 3] += con4;
		}
	}

	delete[] conv1;
	delete[] conv2;

	delete[] filter1o;
	delete[] filter1e;
	delete[] filter2o;
	delete[] filter2e;

	delete[] xe;
	delete[] t;
}


void imageTools::dtwavexfm2(double in[], int H, int W, double out[], double LoLo[], int level)
{
	if (level == 1)
	{
		double* Lo1 = new double[W * H]();
		double* Hi1 = new double[W * H]();
		double* Yh1 = new double[H * W]();

		colfilter(in, H, W, h0o, sh0o, Lo1, false);		//Lo = colfilter(X,h0o).';
		colfilter(in, H, W, h1o, sh1o, Hi1, false);		//Hi = colfilter(X,h1o).';
		if (LoLo)
			colfilter(Lo1, W, H, h0o, sh0o, LoLo, false);		//LoLo = colfilter(Lo, h0o).';

		colfilter(Hi1, W, H, h0o, sh0o, Yh1, false);		//Yh16 = colfilter(Hi,h0o).'
		q2c(Yh1, H, W, out, 1, 6);
		colfilter(Lo1, W, H, h1o, sh1o, Yh1, false);		//Yh34 = colfilter(Lo,h1o).';
		q2c(Yh1, H, W, out, 3, 4);
		colfilter(Hi1, W, H, h1o, sh1o, Yh1, false);		//Yh25 = colfilter(Hi,h1o).';
		q2c(Yh1, H, W, out, 2, 5);

		delete[] Lo1;
		delete[] Hi1;
		delete[] Yh1;
	}
	else if (level >= 2)
	{
		int H2 = H / 2;	//136/2=68
		int W2 = W / 2;	//240/2=120
		double* Lo3 = new double[2 * W2 * H2]();//H=136 W=240
		double* Hi3 = new double[2 * W2 * H2]();
		double* Yh3 = new double[H2 * W2]();

		coldfilt(in, H, W, h0b, h0a, Lo3);		//Lo2 = coldfilt(LoLo1,h0b,h0a).';
		coldfilt(in, H, W, h1b, h1a, Hi3);		//Hi2 = coldfilt(LoLo1,h1b,h1a).';
		coldfilt(Lo3, 2 * W2, H2, h0b, h0a, LoLo);		//LoLo2 = coldfilt(Lo2, h0b, h0a).'

		coldfilt(Hi3, 2 * W2, H2, h0b, h0a, Yh3);		//Yh16 = coldfilt(Hi2,h0b,h0a).'
		q2c(Yh3, H2, W2, out, 1, 6);
		coldfilt(Lo3, 2 * W2, H2, h1b, h1a, Yh3);		//Yh34 = coldfilt(Lo2,h1b,h1a).'
		q2c(Yh3, H2, W2, out, 3, 4);
		coldfilt(Hi3, 2 * W2, H2, h1b, h1a, Yh3);		//Yh25 = coldfilt(Hi2,h1b,h1a).'
		q2c(Yh3, H2, W2, out, 2, 5);

		delete[] Lo3;
		delete[] Hi3;
		delete[] Yh3;
	}
}


void imageTools::dtwaveifm2(double in[], int H, int W, double Zin[], double Zout[], int Level)
{
	if (Level >= 2)
	{
		int H2 = H * 2;//17*2=34
		int W2 = W * 2;//30*2=60
		double *lh3 = new double[H2*W2]();
		double *hl3 = new double[H2*W2]();
		double *hh3 = new double[H2*W2]();
		double *Y13 = new double[W2 * H2 * 2]();	//60*68
		double *Y23 = new double[W2 * H2 * 2]();
		c2q(in, H, W * 2, lh3, 1, 6);		        //lh = c2q(Yh{ current_level }(:, : , [1 6]), gain_mask([1 6], current_level)); = 17, 60,
		c2q(in, H, W * 2, hl3, 3, 4);		        //hl = c2q(Yh{ current_level }(:, : , [3 4]), gain_mask([3 4], current_level));
		c2q(in, H, W * 2, hh3, 2, 5);		        //hh = c2q(Yh{ current_level }(:, : , [2 5]), gain_mask([2 5], current_level));

		colifilt(Zin, H2, W2, g0b, g0a, Y13);	//colifilt(Z, g0b, g0a).' 转置 +=     =34, 60,
		colifilt(lh3, H2, W2, g1b, g1a, Y13);	//colifilt(lh,g1b, g1a).'
		colifilt(hl3, H2, W2, g0b, g0a, Y23);	//colifilt(hl, g0b, g0a).' 转置 +=       =34, 60,
		colifilt(hh3, H2, W2, g1b, g1a, Y23);	//colifilt(hh, g1b, g1a).'     =34, 60,

		memset(Zout, 0, sizeof(double)* H2 * 2 * W2 * 2);         //68 * 120
		colifilt(Y13, W2, H2 * 2, g0b, g0a, Zout);	//colifilt(y1, g0b, g0a).' 转置 +=      = 60, 68,
		colifilt(Y23, W2, H2 * 2, g1b, g1a, Zout);	//colifilt(y2, g1b, g1a).'              = 60, 68,

		delete[] lh3;
		delete[] hl3;
		delete[] hh3;
		delete[] Y13;
		delete[] Y23;
	}
	else if (Level == 1)
	{
		int H2 = H * 2;//68*2=136
		int W2 = W * 2;//120*2=240
		double *lh1 = new double[H2*W2]();
		double *hl1 = new double[H2*W2]();
		double *hh1 = new double[H2*W2]();
		double *Y11 = new double[W2*H2]();	//60*68
		double *Y21 = new double[W2*H2]();
		c2q(in, H, W * 2, lh1, 1, 6);		//lh = c2q(Yh{ current_level }(:, : , [1 6]), gain_mask([1 6], current_level));  //68, 240,
		c2q(in, H, W * 2, hl1, 3, 4);	    //hl = c2q(Yh{ current_level }(:, : , [3 4]), gain_mask([3 4], current_level));
		c2q(in, H, W * 2, hh1, 2, 5);		//hh = c2q(Yh{ current_level }(:, : , [2 5]), gain_mask([2 5], current_level));

		colfilter(Zin, H2, W2, g0o, sg0o, Y11, true); //colfilter(Z, g0o).' = 136, 240,
		colfilter(lh1, H2, W2, g1o, sg1o, Y11, true);	//colfilter(lh, g1o).'
		colfilter(hl1, H2, W2, g0o, sg0o, Y21, true); //colfilter(hl, g0o).'
		colfilter(hh1, H2, W2, g1o, sg1o, Y21, true);	//colfilter(hh, g1o).'

		memset(Zout, 0, sizeof(double)* H2 * W2);       // =  136 * 240
		colfilter(Y11, W2, H2, g0o, sg0o, Zout, true);  //colfilter(y1, g0o).'
		colfilter(Y21, W2, H2, g1o, sg1o, Zout, true);	//colfilter(y2, g1o).'

		delete[] lh1;
		delete[] hl1;
		delete[] hh1;
		delete[] Y11;
		delete[] Y21;

		//以上测试完毕
		/*for (int i = 0; i < 136; i++)
		{
		for (int j = 0; j < 240; j++)
		cout << *(Z1 + i * 240 + j) << " ";
		cout << endl;
		}*/
	}

}


void imageTools::canbe4(vector<vector<double> > arr, vector<vector<double> > *out)
{
	vector<vector<double> > temp4;
	temp4.assign(arr.begin(), arr.end());
	int h = arr.size();
	int w = arr[0].size();

	if (h % 4 != 0)
	{
		extendH(&temp4, 1);
	}
	if (w % 4 != 0)
	{
		extendW(&temp4);
	}
	(*out).assign(temp4.begin(), temp4.end());
	temp4.clear();
	vector<vector<double> >().swap(temp4);

}


void imageTools::extendH(vector<vector<double> > *arr, int H)
{
	double temp_UextendL1;
	vector<double> temp_newUextendL1;
	int arr_size = (*arr).size();
	for (int x = arr_size - 1; x > arr_size - H - 1; x--)         //矩阵最后镜像扩充H行
	{
		for (unsigned int y = 0; y < (*arr)[0].size(); y++)
		{
			temp_UextendL1 = (*arr)[x][y];
			temp_newUextendL1.push_back(temp_UextendL1);
		}
		(*arr).push_back(temp_newUextendL1);
		temp_newUextendL1.clear();
	}


	double temp_UextendL2;
	vector<double> temp_newUextendL2;
	for (int x = 0; x < 2 * H - 1; x = x + 2)                //矩阵最前镜像扩充H行
	{
		for (int y = 0; y < (*arr)[0].size(); y++)
		{
			temp_UextendL2 = (*arr)[x][y];
			temp_newUextendL2.push_back(temp_UextendL2);
		}
		(*arr).insert((*arr).begin(), temp_newUextendL2);
		temp_newUextendL2.clear();
	}


	vector<double>().swap(temp_newUextendL1);
	vector<double>().swap(temp_newUextendL2);
}


void imageTools::extendW(vector<vector<double> > *arr)     //只镜像添加前后一列
{
	double temp_UW;                                   //镜像扩充最后一列
	for (int i = 0; i < (*arr).size(); i++)
	{
		temp_UW = (*arr)[i][(*arr)[i].size() - 1];
		(*arr)[i].push_back(temp_UW);                //相当于在每行最后添加一个数
	}
	double temp_UW1;                                 //镜像扩充最前一列
	for (int x = 0; x < (*arr).size(); x++)
	{
		temp_UW1 = (*arr)[x][0];
		(*arr)[x].insert((*arr)[x].begin(), temp_UW1);
	}
}


void imageTools::getNewRow(vector<vector<double> > in, vector<vector<double> > *out)   //取第2行到 最后一行-1
{
	vector<double> row_temp;
	double temp = 0;
	for (int x = 1; x < in.size() - 1; x++)           //Z = Z(2:row_size-1,:);
	{
		for (int y = 0; y < in[0].size(); y++)
		{
			temp = in[x][y];
			row_temp.push_back(temp);
		}
		(*out).push_back(row_temp);
		row_temp.clear();
	}
	vector<double>().swap(row_temp);
}


void imageTools::getNewCol(vector<vector<double> > in, vector<vector<double> > *out)  //取第二列 到倒数第二列
{
	vector<double> col_temp;
	double temp = 0;
	for (int x = 0; x < in.size(); x++)           //Z = Z(:,2:col_size-1);
	{
		for (int y = 1; y < in[0].size() - 1; y++)
		{
			temp = in[x][y];
			col_temp.push_back(temp);
		}
		(*out).push_back(col_temp);
		col_temp.clear();
	}
	vector<double>().swap(col_temp);
}


void imageTools::getNewRow2(vector<vector<double> > in, vector<vector<double> > *out)   //取第1行到 最后一行-1
{
	vector<double> row_temp;
	double temp = 0;
	for (int x = 0; x < in.size() - 1; x++)           //Z = Z(1:row_size-1,:);
	{
		for (int y = 0; y < in[0].size(); y++)
		{
			temp = in[x][y];
			row_temp.push_back(temp);
		}
		(*out).push_back(row_temp);
		row_temp.clear();
	}
	vector<double>().swap(row_temp);
}


void imageTools::getNewCol2(vector<vector<double> > in, vector<vector<double> > *out)  //取第1列 到倒数第二列
{
	vector<double> col_temp;
	double temp = 0;
	for (int x = 0; x < in.size(); x++)           //Z = Z(:,1:col_size-1);
	{
		for (int y = 0; y < in[0].size() - 1; y++)
		{
			temp = in[x][y];
			col_temp.push_back(temp);
		}
		(*out).push_back(col_temp);
		col_temp.clear();
	}
	vector<double>().swap(col_temp);
}


void imageTools::q2c(double in[], int H, int W, double out[], int  num1, int num2)
{
	int xa, ya, out1, out2;
	double a, b, c, d;
	int H2 = H / 2;
	int W2 = W / 2;
	for (int i = 0; i < H2; i++)
	{
		out1 = (num1 - 1) * H2 * W + i * W;
		out2 = (num2 - 1) * H2 * W + i * W;
		for (int j = 0; j < W2; j++)
		{
			xa = i * 2;
			ya = j * 2;
			a = in[xa * W + ya];
			b = in[xa * W + ya + 1];
			c = in[(xa + 1) * W + ya];
			d = in[(xa + 1) * W + ya + 1];

			out[out1 + j] = (a - d) * sqrt2;   //a-d out1r
			out[out1 + j + W2] = (b + c) * sqrt2;	//b+c out1i
			out[out2 + j] = (a + d) * sqrt2;	//a+d out2r
			out[out2 + j + W2] = (b - c) * sqrt2;	//b-c out2i
		}
	}
}


void imageTools::c2q(double in[], int H, int W, double out[], int  num1, int num2)
{
	int W2 = W / 2;
	int in1, in2;
	for (int i = 0; i < H; i++)
	{
		in1 = (num1 - 1) * H * W + i * W;
		in2 = (num2 - 1) * H * W + i * W;
		for (int j = 0; j < W2; j++)
		{
			out[2 * i * W + 2 * j] = (in[in1 + j] + in[in2 + j]) * sqrt2;
			out[2 * i * W + 2 * j + 1] = (in[in1 + j + W2] + in[in2 + j + W2]) * sqrt2;
			out[(2 * i + 1) * W + 2 * j] = (in[in1 + j + W2] - in[in2 + j + W2]) * sqrt2;
			out[(2 * i + 1) * W + 2 * j + 1] = -(in[in1 + j] - in[in2 + j]) * sqrt2;

		}
	}
}

// 指针已经全部删除
void imageTools::AddWaterMark(BYTE *InU, int *OutU)
{
	double *orgU = new double[m_orgwidth*m_orgheight]();
	for (int i = 0; i < m_orgheight*m_orgwidth; i++)
	{
		orgU[i] = (double)InU[i] - 127;
	}

	int H = m_orgheight;
	int W = m_orgwidth;

	///先看能否整除2，否则添加一行或一列
	vector<vector<double> > intemp;                        //转为二维数组
	double itemp = 0;
	vector<double> temp_in;
	for (int x = 0; x < H; x++)
	{
		for (int y = 0; y < W; y++)
		{
			itemp = orgU[x*W + y];
			temp_in.push_back(itemp);
		}
		intemp.push_back(temp_in);
		temp_in.clear();
	}
	vector<double>().swap(temp_in);

	delete[] orgU;

	if (H % 2 != 0)
	{
		double temp_H;                                     //新添加一行
		vector<double> temp_newUH;
		for (int x = H - 1; x < H; x++)
		{
			for (int y = 0; y < W; y++)
			{
				temp_H = intemp[x][y];
				temp_newUH.push_back(temp_H);
			}
			intemp.push_back(temp_newUH);
			temp_newUH.clear();
		}
		vector<double>().swap(temp_newUH);
	}
	if (W % 2 != 0)                   //新添加一列
	{
		double temp_W;
		for (int i = 0; i < intemp.size(); i++)
		{
			temp_W = intemp[i][W - 1];
			intemp[i].push_back(temp_W);                        //相当于在每行最后添加一个数
		}
	}

	int gai_H = intemp.size();
	int gai_W = intemp[0].size();


	//赋值成员变量
	mgai_H = gai_H;                //1080原始L1L1大小 能否除以2扩充后的 也是colifter之后
	mgai_W = gai_W;                //1920

	double *gai_intemp = new double[gai_H*gai_W]();
	//转回一维数组
	int m = 0;
	for (int x = 0; x < intemp.size(); x++)
	{
		for (int y = 0; y < intemp[0].size(); y++)
		{
			gai_intemp[m] = intemp[x][y];
			m = m + 1;
		}
	}



	L1L1 = new double[gai_H*gai_W]();
	Uh1 = new double[gai_H*gai_W * 3]();
	dtwavexfm2(gai_intemp, gai_H, gai_W, Uh1, L1L1, 1);

	delete[] gai_intemp;


	//判断L1L1的高和宽能否被4整除
	vector<vector<double> > L1L1temp;                        //转为二维数组
	double L1L1itemp = 0;
	vector<double> L1L1temp_in;
	for (int x = 0; x < gai_H; x++)
	{
		for (int y = 0; y < gai_W; y++)
		{
			L1L1itemp = L1L1[x*gai_W + y];
			L1L1temp_in.push_back(L1L1itemp);
		}
		L1L1temp.push_back(L1L1temp_in);
		L1L1temp_in.clear();
	}
	vector<double>().swap(L1L1temp_in);


	vector<vector<double> > tempL1L1;
	canbe4(L1L1temp, &tempL1L1);


	//这时候L1L1这些都是1920*1080(宽*高)  接下来L2 H2为540*1920  L2L2为960*540
	//定义新的宽高 分别为改后宽高的一半
	int newL1L1_H = tempL1L1.size();
	int newL1L1_W = tempL1L1[0].size();
	int gai_Hhalf = newL1L1_H / 2;            //540
	int gai_Whalf = newL1L1_W / 2;            //960

											  //给成员变量赋值 好用于U通道返回时赋值
	mgai_Hhalf = gai_Hhalf;            //540
	mgai_Whalf = gai_Whalf;            //960
	mnewL1L1_H = newL1L1_H;            //1080          //这个是判断L1L1得到后判断能否被4整除后的新L1L1大小
	mnewL1L1_W = newL1L1_W;            //1920

	double *gai_L1L1 = new double[newL1L1_H*newL1L1_W]();
	//转回一维数组
	int m1 = 0;
	for (int x = 0; x < tempL1L1.size(); x++)
	{
		for (int y = 0; y < tempL1L1[0].size(); y++)
		{
			gai_L1L1[m1] = tempL1L1[x][y];
			m1 = m1 + 1;
		}
	}
	vector<vector<double> >().swap(tempL1L1);


	L2L2 = new double[gai_Hhalf*gai_Whalf]();
	Uh2 = new double[gai_Hhalf*gai_Whalf * 3]();

	dtwavexfm2(gai_L1L1, mnewL1L1_H, mnewL1L1_W, Uh2, L2L2, 2);

	delete[] gai_L1L1;


	//判断L2L2的高和宽能否被4整除
	vector<vector<double> > L2L2temp;                        //转为二维数组
	double L2L2itemp = 0;
	vector<double> L2L2temp_in;
	for (int x = 0; x < gai_Hhalf; x++)
	{
		for (int y = 0; y < gai_Whalf; y++)
		{
			L2L2itemp = L2L2[x*gai_Whalf + y];
			L2L2temp_in.push_back(L2L2itemp);
		}
		L2L2temp.push_back(L2L2temp_in);
		L2L2temp_in.clear();
	}
	vector<double>().swap(L2L2temp_in);


	vector<vector<double> > tempL2L2;
	canbe4(L2L2temp, &tempL2L2);


	//这时候L2L2这些都是960*540(宽*高)  接下来L3 H3为270*960  L2L2为480*270
	//定义新的宽高 分别为改后宽高的一半
	int newL2L2_H = tempL2L2.size();          //540
	int newL2L2_W = tempL2L2[0].size();       //960
	int gai_Hhalf2 = newL2L2_H / 2;            //270
	int gai_Whalf2 = newL2L2_W / 2;            //480


											   //赋值成员变量
	mgai_Hhalf2 = gai_Hhalf2;            //270  这是U第三层修改后此时的大小
	mgai_Whalf2 = gai_Whalf2;            //480
	mnewL2L2_H = newL2L2_H;              //540    U第三层往第二层返回时，应有的大小
	mnewL2L2_W = newL2L2_W;              //960



	double *gai_L2L2 = new double[newL2L2_H*newL2L2_W]();
	//转回一维数组
	int m2 = 0;
	for (int x = 0; x < tempL2L2.size(); x++)
	{
		for (int y = 0; y < tempL2L2[0].size(); y++)
		{
			gai_L2L2[m2] = tempL2L2[x][y];
			m2 = m2 + 1;
		}
	}
	vector<vector<double> >().swap(tempL2L2);



	L3L3 = new double[gai_Hhalf2*gai_Whalf2]();
	Uh3 = new double[gai_Hhalf2*gai_Whalf2 * 3]();            //第三层高频m_hblock16 * m_wblock16 * 2 * 6

	dtwavexfm2(gai_L2L2, mnewL2L2_H, mnewL2L2_W, Uh3, L3L3, 3);

	delete[] gai_L2L2;

	int h = gai_Hhalf2 / 2;
	int w = gai_Whalf2 / 2;

	//修改第三层数据的代码
	add2Dtcwt(watermark, h, w, Uh3);


	//恢复
	Z3 = new double[mnewL2L2_W*mnewL2L2_H]();

	dtwaveifm2(Uh3, h, w, L3L3, Z3, 3);


	vector<vector<double> > newL2L2temp;                        //转为二维数组
	double newL2L2itemp = 0;
	vector<double> newL2L2temp_in;
	for (int x = 0; x < mnewL2L2_H; x++)
	{
		for (int y = 0; y < mnewL2L2_W; y++)
		{
			newL2L2itemp = Z3[x*mnewL2L2_W + y];
			newL2L2temp_in.push_back(newL2L2itemp);
		}
		newL2L2temp.push_back(newL2L2temp_in);
		newL2L2temp_in.clear();
	}
	vector<double>().swap(newL2L2temp_in);



	//    int mgai_Hhalf;            //540    原本L2L2的高和宽
	//    int mgai_Whalf;            //960


	if (mgai_Hhalf != mnewL2L2_H)
	{
		vector<vector<double> > newrow;
		getNewRow(newL2L2temp, &newrow);
		newL2L2temp.clear();
		newL2L2temp.assign(newrow.begin(), newrow.end());
		vector<vector<double> >().swap(newrow);
	}
	if (mgai_Whalf != mnewL2L2_W)
	{
		vector<vector<double> > newcol;
		getNewCol(newL2L2temp, &newcol);
		newL2L2temp.clear();
		newL2L2temp.assign(newcol.begin(), newcol.end());
		vector<vector<double> >().swap(newcol);
	}


	//改后newL2L2temp 的宽应为mgai_Whalf  高应为mgai_Hhalf
	double *gai_newL2L2 = new double[mgai_Whalf*mgai_Hhalf]();
	//转回一维数组
	int mL2L2 = 0;
	for (int x = 0; x < newL2L2temp.size(); x++)
	{
		for (int y = 0; y < newL2L2temp[0].size(); y++)
		{
			gai_newL2L2[mL2L2] = newL2L2temp[x][y];
			mL2L2 = mL2L2 + 1;
		}
	}
	vector<vector<double> >().swap(newL2L2temp);

	//    mgai_Hhalf=gai_Hhalf;            //540
	//    mgai_Whalf=gai_Whalf;            //960
	//    mnewL1L1_H=newL1L1_H;            //1080          //这个是判断L1L1得到后判断能否被4整除后的新L1L1大小
	//    mnewL1L1_W=newL1L1_W;            //1920



	Z2 = new double[mnewL1L1_W*mnewL1L1_H]();
	dtwaveifm2(Uh2, mgai_Hhalf / 2, mgai_Whalf / 2, gai_newL2L2, Z2, 2);

	delete[] gai_newL2L2;

	vector<vector<double> > newL1L1temp;                        //转为二维数组
	double newL1L1itemp = 0;
	vector<double> newL1L1temp_in;
	for (int x = 0; x < mnewL1L1_H; x++)
	{
		for (int y = 0; y < mnewL1L1_W; y++)
		{
			newL1L1itemp = Z2[x*mnewL1L1_W + y];
			newL1L1temp_in.push_back(newL1L1itemp);
		}
		newL1L1temp.push_back(newL1L1temp_in);
		newL1L1temp_in.clear();
	}
	vector<double>().swap(newL1L1temp_in);



	//    mgai_H=gai_H;                //1080原始L1L1大小 能否除以2扩充后的 也是colifter之后
	//    mgai_W=gai_W;                //1920

	if (mgai_H != newL1L1temp.size())
	{
		vector<vector<double> > newrow2;
		getNewRow(newL1L1temp, &newrow2);
		newL1L1temp.clear();
		newL1L1temp.assign(newrow2.begin(), newrow2.end());
		vector<vector<double> >().swap(newrow2);
	}
	if (mgai_W != newL1L1temp[0].size())
	{
		vector<vector<double> > newcol2;
		getNewCol(newL1L1temp, &newcol2);
		newL1L1temp.clear();
		newL1L1temp.assign(newcol2.begin(), newcol2.end());
		vector<vector<double> >().swap(newcol2);
	}


	//改后newL1L1temp 的宽应为mgai_W  高应为mgai_H
	double *gai_newL1L1 = new double[mgai_W*mgai_H]();
	//转回一维数组
	int mL1L1 = 0;
	for (int x = 0; x < newL1L1temp.size(); x++)
	{
		for (int y = 0; y < newL1L1temp[0].size(); y++)
		{
			gai_newL1L1[mL1L1] = newL1L1temp[x][y];
			mL1L1 = mL1L1 + 1;
		}
	}
	vector<vector<double> >().swap(newL1L1temp);




	Z1 = new double[mgai_W*mgai_H]();
	dtwaveifm2(Uh1, mgai_H / 2, mgai_W / 2, gai_newL1L1, Z1, 1);

	delete[] gai_newL1L1;

	vector<vector<double> > newU2temp;                        //转为二维数组
	double newU2itemp = 0;
	vector<double> newU2temp_in;
	for (int x = 0; x < mgai_H; x++)
	{
		for (int y = 0; y < mgai_W; y++)
		{
			newU2itemp = Z1[x*mgai_W + y];
			newU2temp_in.push_back(newU2itemp);
		}
		newU2temp.push_back(newU2temp_in);
		newU2temp_in.clear();
	}
	vector<double>().swap(newU2temp_in);


	if (mgai_H != m_orgheight)
	{
		vector<vector<double> > newrow;
		getNewRow2(newU2temp, &newrow);
		newU2temp.clear();
		newU2temp.assign(newrow.begin(), newrow.end());
		vector<vector<double> >().swap(newrow);
	}
	if (m_orgwidth != mgai_W)
	{
		vector<vector<double> > newcol;
		getNewCol2(newU2temp, &newcol);
		newU2temp.clear();
		newU2temp.assign(newcol.begin(), newcol.end());
		vector<vector<double> >().swap(newcol);
	}

	newU2 = new double[m_orgheight*m_orgwidth]();
	int mnewU = 0;
	for (int x = 0; x < m_orgheight; x++)
	{
		for (int y = 0; y < m_orgwidth; y++)
		{
			newU2[mnewU] = newU2temp[x][y];
			mnewU = mnewU + 1;
		}
	}
	vector<vector<double> >().swap(newU2temp);



	double temp = 0;
	for (int x = 0; x < m_orgwidth*m_orgheight; x++)
	{
		temp = newU2[x];
		if (temp > 128)
		{
			newU2[x] = 128;
		}
		if (temp < -127)
		{
			newU2[x] = -127;
		}
		newU2[x] = std::round(newU2[x]) + 127;
		OutU[x] = newU2[x];
	}


}

int imageTools::WaterEncodeDeinit()
{
	delete[] L1L1;
	delete[] L2L2;
	delete[] L3L3;
	delete[] newU2;
	delete[] wh;
	delete[] Uh3;
	delete[] Uhabs;
	delete[] Uh1;
	delete[] Uh2;
	delete[] Z1;
	delete[] Z2;
	delete[] Z3;
	return 0;
}

void imageTools::add2Dtcwt(int w, int H, int W, double Uh1[]){
	vector<Matrix<complex<double>>> Uh = imageTools::rearrange(Uh1, H, W);
	int loc1 = 0;
	int block = 8;
	vector<Matrix<complex<double>>> Um;
	vector<Matrix<complex<double>>> Vm;
	vector<Matrix<double>> Sm;

	//分块后的Sm矩阵
	int Hd, Wd;
	if (H % block == 0)
		Hd = H;
	else
		Hd = (H / block) * block;
	if (W % block == 0)
		Wd = W;
	else
		Wd = (W / block) * block;

	for (int i = 0; i < 6; i++) {
		Matrix< complex<double> > tempA(Hd, Wd);
		Matrix< complex<double> > tempB(Hd, Wd);
		Matrix< double > tempC(Hd, Wd);
		Um.push_back(tempA);
		Vm.push_back(tempB);
		Sm.push_back(tempC);
	}

	//SVD
	for (int i = 0; i < H; i += block) {
		if (i + block > H) {
			continue;
		}
		else {
			for (int j = 0; j < W; j += block) {
				if (j + block > W) {
					continue;
				}
				else {
					for (int k = 0; k < 6; k++) {
						Matrix<complex<double>> A(block, block);
						A = getSubMatrix(Uh[k], i, i + block - 1, j, j + block - 1);
						CSVD<double> svd;
						svd.dec(A);
						Matrix<complex<double>> U = svd.getU();
						Matrix<complex<double>> V = svd.getV();
						Matrix<double> S = svd.getSM();
						Um[k] = setMatrixVal(Um[k], U, i, i + block - 1, j, j + block - 1);
						Vm[k] = setMatrixVal(Vm[k], V, i, i + block - 1, j, j + block - 1);
						Sm[k] = setMatrixVal(Sm[k], S, i, i + block - 1, j, j + block - 1);
					}
					loc1++;
				}
			}
		}
	}

	//取SVD对角阵的第一个系数
	vector<Matrix<double>> DC;
	for (int i = 0; i < 6; i++) {
		Matrix<double> tempDC(Hd / 8, Wd / 8);
		DC.push_back(tempDC);
	}
	loc1 = 0;
	int loc2;
	for (int i = 0; i < H; i += block) {
		if (i + block > H) {
			continue;
		}
		else {
			loc2 = 0;
			for (int j = 0; j < W; j += block) {
				if (j + block > W) {
					continue;
				}
				else {
					for (int k = 0; k < 6; k++) {
						DC[k][loc1][loc2] = Sm[k][i][j];
					}
					loc2++;
				}
			}
			loc1++;
		}
	}

	//取SVD系数的均值
	double* meanDC = new double[6];
	for (int i = 0; i < 6; i++) {
		meanDC[i] = getMatrixMean(DC[i]);
	}

	//开始嵌入
	double delta, temp1, temp2, times1, times2;
	for (int i = 0; i < H; i += block) {
		if (i + block > H) {
			continue;
		}
		else {
			for (int j = 0; j < W; j += block) {
				if (j + block > W) {
					continue;
				}
				else {
					if ((Sm[0][i][j] < 1 / lamda * meanDC[0]) || (Sm[0][i][j] > lamda * meanDC[0]) ||
						(Sm[5][i][j] < 1 / lamda * meanDC[5]) || (Sm[5][i][j] > lamda * meanDC[5]) ||
						(Sm[2][i][j] < 1 / lamda * meanDC[2]) || (Sm[2][i][j] > lamda * meanDC[2]) ||
						(Sm[3][i][j] < 1 / lamda * meanDC[3]) || (Sm[3][i][j] > lamda * meanDC[3]) ||
						(Sm[0][i][j] / Sm[5][i][j] > a) || (Sm[0][i][j] / Sm[5][i][j] < 1 / a) ||
						(Sm[2][i][j] / Sm[3][i][j] > a) || (Sm[2][i][j] / Sm[3][i][j] < 1 / a)) {
						continue;
					}
					else {
						if ((Sm[0][i][j] - Sm[2][i][j] < T) && watermark == 1) {
							delta = T - (Sm[0][i][j] - Sm[2][i][j]);
							temp1 = Sm[0][i][j] + 0.5 * delta;
							temp2 = Sm[2][i][j] - 0.5 * delta;
							times1 = temp1 / Sm[0][i][j];
							times2 = temp2 / Sm[2][i][j];
							Sm[0][i][j] = times1 * Sm[0][i][j];
							Sm[2][i][j] = times2 * Sm[2][i][j];
						}
						else if ((Sm[0][i][j] - Sm[2][i][j] > -T) && watermark == 0) {
							delta = (Sm[0][i][j] - Sm[2][i][j]) + T;
							temp1 = Sm[0][i][j] - 0.5 * delta;
							temp2 = Sm[2][i][j] + 0.5 * delta;
							times1 = temp1 / Sm[0][i][j];
							times2 = temp2 / Sm[2][i][j];
							Sm[0][i][j] = times1 * Sm[0][i][j];
							Sm[2][i][j] = times2 * Sm[2][i][j];
						}
						if ((Sm[5][i][j] - Sm[3][i][j] < T) && watermark == 1) {
							delta = T - (Sm[5][i][j] - Sm[3][i][j]);
							temp1 = Sm[5][i][j] + 0.5 * delta;
							temp2 = Sm[3][i][j] - 0.5 * delta;
							times1 = temp1 / Sm[5][i][j];
							times2 = temp2 / Sm[3][i][j];
							Sm[5][i][j] = times1 * Sm[5][i][j];
							Sm[3][i][j] = times2 * Sm[3][i][j];
						}
						else if ((Sm[5][i][j] - Sm[3][i][j] > -T) && watermark == 0) {
							delta = (Sm[5][i][j] - Sm[3][i][j]) + T;
							temp1 = Sm[5][i][j] - 0.5 * delta;
							temp2 = Sm[3][i][j] + 0.5 * delta;
							times1 = temp1 / Sm[5][i][j];
							times2 = temp2 / Sm[3][i][j];
							Sm[5][i][j] = times1 * Sm[5][i][j];
							Sm[3][i][j] = times2 * Sm[3][i][j];
						}
					}
				}
			}
		}
	}
	
	//逆SVD
	vector<Matrix<complex<double>>> Uh3;
	Matrix<complex<double>> tempUm (block, block);
	Matrix<complex<double>> tempVm(block, block);
	Matrix<double> tempSm(block, block);
	Matrix<complex<double>> tempOut(block, block);
	for (int i = 0; i < 6; i++) {
		Matrix<complex<double>> tempUh3(H, W);
		Uh3.push_back(tempUh3);
	}
	for (int i = 0; i < H; i += block) {
		if (i + block > H) {
			continue;
		}
		else {
			for (int j = 0; j < W; j += block) {
				if (j + block > W) {
					continue;
				}
				else {
					for (int k = 0; k < 6; k++) {
						Matrix<complex<double>> tempRes(block, block);
						tempSm = getSubMatrix(Sm[k], i, i + block - 1, j, j + block - 1);
						tempUm = getSubMatrix(Um[k], i, i + block - 1, j, j + block - 1);
						tempVm = getSubMatrix(Vm[k], i, i + block - 1, j, j + block - 1);
						Matrix<complex<double>> CS = complexMatrix(tempSm);
						tempOut = tempUm * CS * trH(tempVm);
						Uh3[k] = setMatrixVal(Uh3[k], tempOut, i, i + block - 1, j, j + block - 1);
					}
				}
			}
		}
	}
	Uh1 = arrange(Uh3, H, W);
	
	delete meanDC;
}

vector<Matrix<complex<double>>> imageTools::rearrange(double Uh[], int H, int W) {
	int H2 = H / 2;
	int W2 = W / 2;
	int out1, out2;
	vector<Matrix<complex<double>>> res;
	for (int i = 0; i < 6; i++) {
		Matrix< complex<double> > A(H, W);
		res.push_back(A);
	}
	for (int num = 0; num < 6; num++) {
		for (int i = 0; i < H; i++) {
			out1 = num * H2 * W + i * W;
			for (int j = 0; j < W; j++) {

				double real = Uh[out1 + j];
				double imag = Uh[out1 + j + W2];
				res[num][i][j] = complex<double>(real, imag);
			}
		}
	}
	return res;
}


double* imageTools::arrange(vector<Matrix<complex<double>>> res, int H, int W) {
	double* out = new double[res.size()*res[0].cols()*res[0].rows() * 2];

	int H2 = H / 2;
	int W2 = W / 2;
	int out1;
	for (int num = 0; num < 6; num++) {
		for (int i = 0; i < H2; i++)
		{
			out1 = num * H2 * W + i * W;
			for (int j = 0; j < W2; j++)
			{
				out[out1 + j] = res[num][i][j].real();   //a-d out1r
				out[out1 + j + W2] = res[num][i][j].imag(); //b+c out1i
			}
		}
	}

	return out;
}

double imageTools::getMatrixMean(Matrix<double> in) {
	double sum = 0;
	for (int i = 0; i < in.rows(); i++) {
		for (int j = 0; j < in.cols(); j++) {
			sum += in[i][j];
		}
	}
	double num = (double)in.rows()*in.cols();
	return sum / num;
}
