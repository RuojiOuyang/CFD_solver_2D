#pragma once

#define PI 3.1415926535

class CCell//CCell类表示中心点
{
public:
	double x;//中心点x坐标
	double y;//中心点y坐标
	double Den;//密度
	double u;//x向速度
	double v;//y向速度
	double E;//总能量
	double H;//总焓
	double P;//压力
	double T;//温度
	double Ma;//马赫数
	double Vol;//控制体体积
	double localdt;//当地时间步长
	double dt;//最小时间步长
	double c;//当地声速
	double U[6];//中心守恒变量，取6是为了可以从1开始，后面一样
	double QC[6];//控制体的对流通量积分（净无粘通量积分）
	double QV[6];//控制体的扩散通量积分（净粘性通量积分）
	double D[6];//控制体的人工粘性通量积分
	double Med[6];//RK法计算过程中的残量（式2-1-12中的R）
	double U_former[6];//用于保存Runge-Kutta中上一次迭代的值
	double T_x;//T在x方向的偏导数
	double T_y;//T在y方向的偏导数
	double u_x;//速度u在x方向的偏导数
	double u_y;//速度u在y方向的偏导数
	double v_x;//速度v在x方向的偏导数
	double v_y;//速度v在y方向的偏导数
	double Touxx;//剪切力Touxx项
	double Touxy;//剪切力Touxy项
	double Touyy;//剪切力Touyy项
	double Touhx;//式（3-1-1）中的Touhx
	double Touhy;//式（3-1-1）中的Touhy
	double miu;//气体动力粘度
	double miut;//涡粘粘度，乘了fv1之后的，即Den*miubl*fv1
	double miubl;//湍流粘度，SA模型中的长得像小v的
	double miubl_x;//湍流粘度在x方向的偏导数
	double miubl_y;//湍流粘度在y方向的偏导数
	double Touvx;//湍流模型Touvx项
	double Touvy;//湍流模型Touvy项
	double S;//源项
	double writeX;//手写体X，即湍流粘度与层流粘度之比
	double fv1;//SA模型系数
	double d;//中心点距最近的物面的距离，SA模型
};

class CFace//CFace类表示控制体界面
{
public:
	double nx;//边的法向面积的x分量，算当地时间步的时候用的
	double ny;//边的法向面积的y分量
	double Den;//边的密度
	double u;//x向速度
	double v;//y向速度
	double E;//边的总能量
	double H;//边的总焓
	double P;//边的压力
	double T;//边的温度
	double FU[6];//FU是各边对流通量变量，既然定义出了这个，就说明我是按照先平均守恒量，再计算通量的方式来计算QC的，而不是先计算通量再平均
	double FF[6];//各边的对流通量与面积矢量的点积F.*S
	double FK[6];//各边的扩散通量
	double FD[6];//各边的人工粘性通量
	double lamta;//粘性项里的系数，计算人工粘性的时候用
	double epslion[3];//粘性项里自适应系数
	double dU[6];//守恒量的一阶差分
	double dddU[6];//守恒量的三阶差分
	double miubl;//湍流粘度，SA模型中的
	double miubl_x;//湍流粘度的偏导数，SA模型中的
	double miubl_y;
};

class CNode//CNode类是装网格点用的
{
public:
	double x;//网格点x坐标
	double y;//网格点y坐标
	double Den;//网格点处密度
	double u;//网格点处x向速度
	double v;//网格点处y向速度
	double E;//网格点处总能量
	double H;//网格点处总焓
	double P;//网格点处压力
	double T;//网格点处温度
	double Ma;//网格点处马赫数
	double v_x;
	double u_y;//这两个量用来计算涡量
};

extern CNode** Node;//网格点信息
extern CCell** Cell;//中心点
extern CCell** DummyIS;//虚元：I方向底端，south
extern CCell** DummyIN;//虚元：I方向顶端，north
extern CCell** DummyJW;//虚元：J方向左端，west
extern CCell** DummyJE;//虚元：J方向右端，east
extern CFace** FaceI;//I向边界，指顺着i方向的边界
extern CFace** FaceJ;//J向边界，指顺着j方向的边界

extern int i, j;//最常见的循环变量
extern int Imax, Jmax;//网格点数，i是周向，j是法向
extern int n_dummy;//虚元（ghost cell）的数量
extern int n;//表示迭代循环的
extern double k;//比热比
extern double R;//气体常数
extern double T0;//海平面温度参考值
extern double Ts;//Sutherland常数
extern double mu0;//对应于T0的动力粘度
extern double P0;//海平面压力参考值
extern double c0;//海平面声速参考值
extern double cv;
extern double cp;
extern double Pr;//普朗特数
extern double CFL;//CFL
extern double AOA;//攻角，单位度
extern int Iternum_max;//最大迭代代数
extern double* RKalpha;//Runge_Kutta的系数
extern double totaltime;//仿真推进的总时间
extern double** Den_former;
extern double** muI;//I方向边界的激波感受因子
extern double** muJ;//J方向边界的激波感受因子
extern double Resid;//残差
extern double Tt, Pt, Dent, ut, vt, Ma, ct, Et, Ht, miut, miublt;//初始来流数据
extern double Cb1, Cb2, delta, writek, rlim, Cw1, Cw2, Cw3, Cv1, Prt;//SA湍流模型中的常数

extern void Dynamic_array();