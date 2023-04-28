/*
变量定义与声明都在这里
类的定义放在同名头文件里
*/

/*
变量定义与声明都在这里
类的定义放在同名头文件里
*/

#include "Main.h"
#include "VariablesDeclear.h"

//1.中心点
arma::mat Cell_Den;//密度
arma::mat Cell_u;//x向速度
arma::mat Cell_v;//y向速度
arma::mat Cell_E;//总能量
arma::mat Cell_P;//压力
arma::mat Cell_T;//温度
arma::mat Cell_Ma;//马赫数
arma::mat Cell_dt;//时间步长
arma::mat Cell_c;//当地声速
arma::mat Cell_U1;//守恒量第一项
arma::mat Cell_U2;//守恒量第二项
arma::mat Cell_U3;//守恒量第三项
arma::mat Cell_U4;//守恒量第四项
arma::mat Cell_U1_former;//前一步守恒量第一项
arma::mat Cell_U2_former;//前一步守恒量第二项
arma::mat Cell_U3_former;//前一步守恒量第三项
arma::mat Cell_U4_former;//前一步守恒量第四项
arma::mat Cell_QC1;//控制体的对流通量积分第一项
arma::mat Cell_QC2;//控制体的对流通量积分第二项
arma::mat Cell_QC3;//控制体的对流通量积分第三项
arma::mat Cell_QC4;//控制体的对流通量积分第四项
arma::mat Cell_D1;//控制体的人工粘性通量积分第一项
arma::mat Cell_D2;//控制体的人工粘性通量积分第二项
arma::mat Cell_D3;//控制体的人工粘性通量积分第三项
arma::mat Cell_D4;//控制体的人工粘性通量积分第四项
arma::mat Cell_Med1;//RK法计算过程中的残量第一项
arma::mat Cell_Med2;//RK法计算过程中的残量第二项
arma::mat Cell_Med3;//RK法计算过程中的残量第三项
arma::mat Cell_Med4;//RK法计算过程中的残量第四项
arma::mat Cell_x;//中心点x坐标
arma::mat Cell_y;//中心点y坐标
arma::mat Cell_Vol;//控制体体积

//2.控制体界面
//2.1I向界面（横着的）
arma::mat FaceI_Den;//密度
arma::mat FaceI_u;//x向速度
arma::mat FaceI_v;//y向速度
arma::mat FaceI_E;//总能量
arma::mat FaceI_P;//压力
arma::mat FaceI_T;//温度
arma::mat FaceI_U1;//界面守恒量第一项
arma::mat FaceI_U2;//界面守恒量第二项
arma::mat FaceI_U3;//界面守恒量第三项
arma::mat FaceI_U4;//界面守恒量第四项
arma::mat FaceI_F1;//各边的对流通量与面积矢量的点积F*S第一项
arma::mat FaceI_F2;//各边的对流通量与面积矢量的点积F*S第二项
arma::mat FaceI_F3;//各边的对流通量与面积矢量的点积F*S第三项
arma::mat FaceI_F4;//各边的对流通量与面积矢量的点积F*S第四项
arma::mat FaceI_lamta;////粘性项里的系数，计算人工粘性的时候用
arma::mat FaceI_epslion1;//粘性项里自适应系数
arma::mat FaceI_epslion2;//粘性项里自适应系数
arma::mat FaceI_dU1;//守恒量的一阶差分第一项
arma::mat FaceI_dU2;//守恒量的一阶差分第二项
arma::mat FaceI_dU3;//守恒量的一阶差分第三项
arma::mat FaceI_dU4;//守恒量的一阶差分第四项
arma::mat FaceI_dddU1;//守恒量的三阶差分第一项
arma::mat FaceI_dddU2;//守恒量的三阶差分第二项
arma::mat FaceI_dddU3;//守恒量的三阶差分第三项
arma::mat FaceI_dddU4;//守恒量的三阶差分第四项
arma::mat FaceI_D1;//各边的人工粘性通量第一项
arma::mat FaceI_D2;//各边的人工粘性通量第二项
arma::mat FaceI_D3;//各边的人工粘性通量第三项
arma::mat FaceI_D4;//各边的人工粘性通量第四项
arma::mat FaceI_nx;//边的法向面积的x分量
arma::mat FaceI_ny;//边的法向面积的y分量
//2.2J向界面（竖着的）
arma::mat FaceJ_Den;//密度
arma::mat FaceJ_u;//x向速度
arma::mat FaceJ_v;//y向速度
arma::mat FaceJ_E;//总能量
arma::mat FaceJ_P;//压力
arma::mat FaceJ_T;//温度
arma::mat FaceJ_U1;//界面守恒量第一项
arma::mat FaceJ_U2;//界面守恒量第二项
arma::mat FaceJ_U3;//界面守恒量第三项
arma::mat FaceJ_U4;//界面守恒量第四项
arma::mat FaceJ_F1;//各边的对流通量与面积矢量的点积F*S第一项
arma::mat FaceJ_F2;//各边的对流通量与面积矢量的点积F*S第二项
arma::mat FaceJ_F3;//各边的对流通量与面积矢量的点积F*S第三项
arma::mat FaceJ_F4;//各边的对流通量与面积矢量的点积F*S第四项
arma::mat FaceJ_lamta;//粘性项里的系数，计算人工粘性的时候用
arma::mat FaceJ_epslion1;//粘性项里自适应系数
arma::mat FaceJ_epslion2;//粘性项里自适应系数
arma::mat FaceJ_dU1;//守恒量的一阶差分第一项
arma::mat FaceJ_dU2;//守恒量的一阶差分第二项
arma::mat FaceJ_dU3;//守恒量的一阶差分第三项
arma::mat FaceJ_dU4;//守恒量的一阶差分第四项
arma::mat FaceJ_dddU1;//守恒量的三阶差分第一项
arma::mat FaceJ_dddU2;//守恒量的三阶差分第二项
arma::mat FaceJ_dddU3;//守恒量的三阶差分第三项
arma::mat FaceJ_dddU4;//守恒量的三阶差分第四项
arma::mat FaceJ_D1;//各边的人工粘性通量第一项
arma::mat FaceJ_D2;//各边的人工粘性通量第二项
arma::mat FaceJ_D3;//各边的人工粘性通量第三项
arma::mat FaceJ_D4;//各边的人工粘性通量第四项
arma::mat FaceJ_nx;//边的法向面积的x分量
arma::mat FaceJ_ny;//边的法向面积的x分量

//3.网格点
arma::mat Node_Den;//密度
arma::mat Node_u;//x向速度
arma::mat Node_v;//y向速度
arma::mat Node_E;//总能量
arma::mat Node_P;//压力
arma::mat Node_T;//温度
arma::mat Node_Ma;//马赫数
arma::mat Node_x;
arma::mat Node_y;

//4.虚元
//4.1I方向底端，south
arma::mat DummyIS_Den;//密度
arma::mat DummyIS_u;//x向速度
arma::mat DummyIS_v;//y向速度
arma::mat DummyIS_E;//总能量
arma::mat DummyIS_P;//压力
arma::mat DummyIS_T;//温度
arma::mat DummyIS_Ma;//马赫数
arma::mat DummyIS_dt;//时间步长
arma::mat DummyIS_c;//当地声速
arma::mat DummyIS_U1;//守恒量第一项
arma::mat DummyIS_U2;//守恒量第二项
arma::mat DummyIS_U3;//守恒量第三项
arma::mat DummyIS_U4;//守恒量第四项
//4.2I方向顶端，north
arma::mat DummyIN_Den;//密度
arma::mat DummyIN_u;//x向速度
arma::mat DummyIN_v;//y向速度
arma::mat DummyIN_E;//总能量
arma::mat DummyIN_P;//压力
arma::mat DummyIN_T;//温度
arma::mat DummyIN_Ma;//马赫数
arma::mat DummyIN_dt;//时间步长
arma::mat DummyIN_c;//当地声速
arma::mat DummyIN_U1;//守恒量第一项
arma::mat DummyIN_U2;//守恒量第二项
arma::mat DummyIN_U3;//守恒量第三项
arma::mat DummyIN_U4;//守恒量第四项
//4.3J方向左端，west
arma::mat DummyJW_Den;//密度
arma::mat DummyJW_u;//x向速度
arma::mat DummyJW_v;//y向速度
arma::mat DummyJW_E;//总能量
arma::mat DummyJW_P;//压力
arma::mat DummyJW_T;//温度
arma::mat DummyJW_Ma;//马赫数
arma::mat DummyJW_dt;//时间步长
arma::mat DummyJW_c;//当地声速
arma::mat DummyJW_U1;//守恒量第一项
arma::mat DummyJW_U2;//守恒量第二项
arma::mat DummyJW_U3;//守恒量第三项
arma::mat DummyJW_U4;//守恒量第四项
//4.4J方向右端，east
arma::mat DummyJE_Den;//密度
arma::mat DummyJE_u;//x向速度
arma::mat DummyJE_v;//y向速度
arma::mat DummyJE_E;//总能量
arma::mat DummyJE_P;//压力
arma::mat DummyJE_T;//温度
arma::mat DummyJE_Ma;//马赫数
arma::mat DummyJE_dt;//时间步长
arma::mat DummyJE_c;//当地声速
arma::mat DummyJE_U1;//守恒量第一项
arma::mat DummyJE_U2;//守恒量第二项
arma::mat DummyJE_U3;//守恒量第三项
arma::mat DummyJE_U4;//守恒量第四项

//网格信息变量
int Imax, Jmax;//网格点数，i是周向，j是法向

//流场信息变量
int n_dummy;//虚元（ghost cell）的数量

//计算过程所需要的变量
int i, j;//最常见的循环变量
double CFL;//CFL数
int n;//当前迭代数
int Iternum_max;//最大迭代代数
double* RKalpha;//Runge_Kutta的系数
double Resid;//残差
arma::mat P_former;//Den_former[]用于保存上一步的各点密度,用于残差计算
arma::mat muI;//I方向的激波感受因子，人工粘性用
arma::mat muJ;//J方向的激波感受因子，人工粘性用
int numAOA;//攻角计数
int M;//配点数

//关于气体的常数
double k;//比热比
double R;//气体常数
double T0;//海平面温度参考值
double P0;//海平面压力参考值
double c0;//海平面声速参考值

//来流数据
double Tt;//来流温度
double Pt;//来流压力
double Dent;//来流密度
double ut;//来流水平速度
double vt;//来流竖直速度
double Ma;//来流马赫数
double ct;//来流声速
double Et;//来流内能
double AOA;//来流攻角，单位度