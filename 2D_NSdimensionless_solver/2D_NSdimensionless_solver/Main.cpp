#include "Main.h"
#include "Read_mesh.h"
#include "Geometry.h"
#include "Initialization.h"
#include "Solve_conservation.h"
#include "Runge_Kutta.h"
#include "Residual.h"
#include "Result_deal.h"
#include "Putout.h"

using namespace std;

int main(int argc, char* argv[])
{
	totaltime = 0;

	ifstream read;
	read.open("yuanzhudata.txt");
	read >> Imax >> Jmax;//读前两个数，前两个数分别是i,j方向网格点的数量
	read.close();
	n_dummy = 3;//虚元的数量
	Dynamic_array();//分配内存

	//初始条件部分
	//常数设定
	k = 1.4;//比热比
	R = 287.06;//气体常数
	T0 = 288.16;//温度参考值
	Ts = 110;//Sutherland常数
	P0 = 101325.0;//压力参考值
	mu0 = 1.7894e-5;//对应于T0的动力粘度
	Pr = 0.71;//普朗特数
	CFL = 2;//CFL
	AOA = 0;//攻角，单位度
	RKalpha[1] = 1.0 / 4.0;//RK法的系数
	RKalpha[2] = 1.0 / 6.0;
	RKalpha[3] = 3.0 / 8.0;
	RKalpha[4] = 0.5;
	RKalpha[5] = 1.0;
	Iternum_max = 10000;//最大迭代步数

	//进口条件
	//L = 50 * 1e-3;//特征尺寸
	Re = 228;//来流雷诺数
	Ma = 0.0002;
	Ttinf = T0 / (1.0 + (k - 1.0) / 2.0 * Ma * Ma);//来流温度
	Ttinf = T0;
	mutinf = mu0 * pow(Ttinf / T0, 1.5) * (T0 + Ts) / (Ttinf + Ts);//来流动力黏度
	//P0 = Re * R * Ttinf * mutinf / (pow(Ttinf / T0, k / (k - 1.0)) * sqrt(k * R * Ttinf) * Ma * L);//求出压力参考值
	Ptinf = P0 * pow(Ttinf / T0, k / (k - 1.0));//来流压力
	ctinf = sqrt(k * R * Ttinf);//来流声速
	Dentinf = Ptinf / (R * Ttinf);//来流密度
	utinf = ctinf * Ma * cos(AOA / 180 * PI);//来流x方向速度
	vtinf = ctinf * Ma * sin(AOA / 180 * PI);//来流y方向速度
	Etinf = Ptinf / (k - 1.0) / Dentinf + (utinf * utinf + vtinf * vtinf) / 2.0;//来流总能量
	Htinf = Etinf + Ptinf / Dentinf;//来流总焓
	L = Re * mutinf / Dentinf / Ma / ctinf;//特征尺寸

	Tt = 1;
	Pt = Ptinf / (Dentinf * ctinf * ctinf);
	ct = 1;
	Dent = 1;
	ut = utinf / ctinf;
	vt = vtinf / ctinf;
	Et = Etinf / ctinf / ctinf;
	Ht = Htinf / ctinf / ctinf;
	mut = 1;

	//网格处理部分
	Read_mesh();//读网格
	Geometry();//处理网格，得到控制体体积以及面积矢量

	//计算部分
	Initialization();//利用上面的初始条件初始化流场

	Solve_conservation();//求各个中心点上的守恒量

	ofstream outx;
	outx.open("Res.dat");//写结果
	for (n = 1; n <= Iternum_max; n++)//求解的循环过程
	{
		Runge_Kutta();//RK法向前求解
		Residual();//计算残差
		cout << n << "   " << Resid << endl;
		outx << n << "   " << Resid << endl;
		if (Resid <= 1e-18)
		{
			break;
		}
	}
	outx.close();

	Result_deal();//处理数据

	Putout();//输出数据

	cout << totaltime * L / ctinf << endl;

	for (n = Iternum_max + 1;; n++)
	{
		Runge_Kutta();//RK法向前求解
		Residual();//计算残差
		cout << n << "   " << Resid << endl;
		outx << n << "   " << Resid << endl;
		if (n % 10000 == 0)
		{
			Result_deal();//处理数据
			Putout();//输出数据
			cout << totaltime * L / ctinf << endl;
		}
	}

	return 0;
}