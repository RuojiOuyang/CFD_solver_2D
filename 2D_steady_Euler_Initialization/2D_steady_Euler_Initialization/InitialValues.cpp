/*
InitialValues函数给出给出初始状态
*/

#include "Main.h"
#include "DynamicArray.h"
#include "ReadMesh.h"
#include "Geometry.h"

using namespace std;
using namespace arma;

void InitialValues()
{
	//1.有关空气的常数
	k = 1.4;//比热比
	R = 287.06;//气体常数
	T0 = 288.16;//海平面温度参考值
	P0 = 101325.0;//海平面压力参考值
	c0 = 340.28;//海平面声速参考值

	//2.来流数据
	double ti;
	double Omega;
	Ma = 0.6;//来流马赫数
	Tt = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//来流温度
	Pt = P0 * pow(Tt / T0, k / (k - 1.0));//来流压力
	ct = sqrt(k * R * Tt);//来流声速
	Dent = Pt / (R * Tt);//来流密度

	//Omega = 2 * ct * Ma * 0.0808 / 1.0;
	//ti = 2 * PI / Omega * ((long long)numAOA - 1.0) / M;//时域配点
	//AOA = 2.89 + 2.41 * sin(Omega * ti);//攻角，单位度
	AOA = 3;

	ut = ct * Ma * cos(AOA / 180 * PI);//来流x方向速度
	vt = ct * Ma * sin(AOA / 180 * PI);//来流y方向速度
	Et = Pt / (k - 1.0) / Dent + (ut * ut + vt * vt) / 2.0;//来流总能量

	//3.接下来需要分配内存，分配内存之前需要一些关键参数：Imax,Jmax,n_dummy
	ifstream read;
	read.open("naca0012_09_17data.txt");
	read >> Imax >> Jmax;//读前两个数，前两个数分别是i,j方向网格点的数量
	read.close();
	n_dummy = 3;
	DynamicArray();

	//计算过程需要的量
	CFL = 2;//CFL数
	RKalpha[1] = 1.0 / 4.0;//RK法的系数
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 10000;//最大迭代步数

	//网格处理部分
	ReadMesh();//读网格
	Geometry();//处理网格，得到控制体体积以及面积矢量
}