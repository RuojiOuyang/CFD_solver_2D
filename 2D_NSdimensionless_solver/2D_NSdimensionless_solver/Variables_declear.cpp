/*
变量定义与声明都在这里
类的定义放在同名头文件里
Dynamic_array函数来预分配内存
*/

#include "Variables_declear.h"

using namespace std;

CNode** Node;//网格点信息
CCell** Cell;//中心点
CCell** DummyIS;//虚元：I方向底端，south
CCell** DummyIN;//虚元：I方向顶端，north
CCell** DummyJW;//虚元：J方向左端，west
CCell** DummyJE;//虚元：J方向右端，east
CFace** FaceI;//I向边界，指顺着i方向的边界
CFace** FaceJ;//J向边界，指顺着j方向的边界

int i, j;//最常见的循环变量
int Imax, Jmax;//网格点数，i是周向，j是法向
int n_dummy;//虚元（ghost cell）的数量
int n;//表示迭代循环的

//关于气体的常数
double k;//比热比
double R;//气体常数
double T0;//海平面温度参考值
double Ts;//Sutherland常数
double mu0;//对应于T0的动力粘度
double P0;//海平面压力参考值
double c0;//海平面声速参考值
double cv;
double cp;
double Pr;//普朗特数

//计算过程的常数
double CFL;//CFL
double AOA;//攻角，单位度
int Iternum_max;//最大迭代代数
double* RKalpha;//Runge_Kutta的系数
double totaltime;//仿真推进的总时间

//计算过程中的变量
double** Den_former;//Den_former[]用于保存上一步的各点密度,用于残差计算
double** muI;//I方向的激波感受因子
double** muJ;//J方向的激波感受因子
double Resid;//残差

//来流数据
double Re;//雷诺数
double Ttinf, Ptinf, Dentinf, utinf, vtinf, Ma, ctinf, Etinf, Htinf, mutinf;//初始来流数据
double Tt, Pt, Dent, ut, vt, ct, Et, Ht, mut;//初始来流数据
double L;//特征尺寸

void Dynamic_array()
{
	//后面写的维数都是列*行，为了与网格的i方向与j方向保持一致
	Node = new CNode * [Imax + 1];//Node是Imax*Jmax维的。+1是为了让标号可以从1开始，后面一样
	Cell = new CCell * [Imax];//Cell是(Imax-1)*(Jmax-1)维的。中心点数量比网格点少1
	DummyIS = new CCell * [Imax];//IS是(Imax-1)*n_dummy维的
	DummyIN = new CCell * [Imax];//IN是(Imax-1)*n_dummy维的
	DummyJW = new CCell * [n_dummy + 1];//JW是n_dummy*(Jmax-1)维的
	DummyJE = new CCell * [n_dummy + 1];//JE是n_dummy*(Jmax-1)维的
	FaceI = new CFace * [Imax];//FaceI是(Imax-1)*Jmax维的边的数量和网格点的数量不一致，实际上，假如网格4*4，则I方向的边是3*4的，J方向是4*3的，总共24条边
	FaceJ = new CFace * [Imax + 1];//FaceJ是Imax*(Jmax-1)维的。边的数量和网格点的数量不一致
	Den_former = new double* [Imax];//Den_former与Cell维数一致，用于保存上一步的各点密度,用于残差计算
	muI = new double* [Imax];//激波感受因子的数量分别和中心点的数量+2*虚元数量一致，竖着算的记为I方向，横着算的记为J方向，muI是(Imax-1)*(Jmax+2*n_dummy)维
	muJ = new double* [Imax + 2 * n_dummy + 1];//muJ是(Imax+2*n_dummy)*(Jmax-1)维的

	for (i = 0; i <= Imax; i++)
	{
		Node[i] = new CNode[Jmax + 1];
		FaceJ[i] = new CFace[Jmax];
	}

	for (i = 0; i <= Imax - 1; i++)
	{
		Cell[i] = new CCell[Jmax];
		FaceI[i] = new CFace[Jmax + 1];
		DummyIS[i] = new CCell[n_dummy + 1];
		DummyIN[i] = new CCell[n_dummy + 1];
		Den_former[i] = new double[Jmax];
		muI[i] = new double[Jmax + 2 * n_dummy + 1];
	}

	for (i = 0; i <= Imax + 2 * n_dummy; i++)
	{
		muJ[i] = new double[Jmax];
	}

	for (i = 0; i <= n_dummy; i++)
	{
		DummyJW[i] = new CCell[Jmax];
		DummyJE[i] = new CCell[Jmax];
	}

	RKalpha = new double[6];//RK的系数
}